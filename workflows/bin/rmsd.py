#!/usr/bin/env python
from tcrtrifold.tcr import annotate_tcr
from mdaf3.AF3OutputParser import AF3Output
from mdaf3.BoltzOutputParser import BoltzOutput
from mdaf3.FeatureExtraction import split_apply_combine
import MDAnalysis as mda
import numpy as np
import argparse
import polars as pl
from pathlib import Path
from MDAnalysis.analysis import rms
import polars.selectors as cs
import Bio.pairwise2


def seq_align_residue_groups(
    mobile_res,
    reference_res,
    mismatch_penalty=-float("inf"),
    match_score=2,
    gap_penalty=-2,
    gapextension_penalty=-0.1,
):
    """
    Alignment = namedtuple("Alignment", ("seqA, seqB, score, start, end"))
    """

    aln = Bio.pairwise2.align.globalms(
        reference_res.sequence(format="string"),
        mobile_res.sequence(format="string"),
        match_score,
        mismatch_penalty,
        gap_penalty,
        gapextension_penalty,
    )[0]

    mobile_aln = aln[1]
    reference_aln = aln[0]

    mobile_res_boolmask = np.full(len(mobile_res), False, dtype=np.bool)
    reference_res_boolmask = np.full(len(reference_res), False, dtype=np.bool)

    true_ptr = 0
    pred_ptr = 0
    for i in range(len(reference_aln)):
        if mobile_aln[i] == "-":
            pred_ptr += 1
            continue
        elif reference_aln[i] == "-":
            true_ptr += 1
            continue
        elif mobile_aln[i] != reference_aln[i]:
            raise ValueError("Sequences cannot align without mutations")
        else:
            mobile_res_boolmask[true_ptr] = True
            reference_res_boolmask[pred_ptr] = True
            true_ptr += 1
            pred_ptr += 1

    # rg1_subseq = rg1[[True if aln_symbol != "-" else False for aln_symbol in rg1_aln]]
    # rg2_subseq = rg2[[True if aln_symbol != "-" else False for aln_symbol in rg2_aln]]

    mobile_res_subseq = mobile_res[mobile_res_boolmask]
    reference_res_subseq = reference_res[reference_res_boolmask]

    assert mobile_res_subseq.sequence(format="string") == reference_res_subseq.sequence(
        format="string"
    )

    return mobile_res_subseq, reference_res_subseq


def align_mhc_to_origin(u, mhc_axes, mhc_origin):

    u_rot = np.linalg.solve(mhc_axes.T, np.eye(3, 3)).T

    # translate and rotate positions in a way that puts
    # mhc coordinate frame at the origin and aligns it
    # to the coordinate axes

    u.atoms.positions = (u.atoms.positions + (mhc_origin * -1)) @ u_rot

    return u


def true_pred_docking_rmsd(row, inf_type, clean_pdb_dir, inf_dir, rank):
    """
    Implementation of CDR center-of-mass RMSD after MHC alignment
    """

    sample_id = row[f"pred_sample_rank_{rank}"]

    if inf_type == "af3":
        inf_cls = AF3Output
        sample_kwarg = {"seed": sample_id}
    elif inf_type == "boltz":
        inf_cls = BoltzOutput
        sample_kwarg = {"sample_num": sample_id}
    else:
        raise ValueError

    pred_u = inf_cls(inf_dir / row["job_name"]).get_mda_universe(**sample_kwarg)
    pred_mhc_axes = np.array(row[f"pred_mhc_axes_{rank}"])
    pred_mhc_origin = np.array(row[f"pred_mhc_origin_{rank}"])
    pred_tcr_axes = np.array(row[f"pred_tcr_axes_{rank}"])
    pred_tcr_origin = np.array(row[f"pred_tcr_origin_{rank}"])
    pred_dgeom = np.array(row[f"pred_dgeom_{rank}"])

    pred_u = align_mhc_to_origin(pred_u, pred_mhc_axes, pred_mhc_origin)

    true_u = mda.Universe(clean_pdb_dir / f"{row['pdb']}.pdb")
    true_mhc_axes = np.array(row[f"true_mhc_axes"])
    true_mhc_origin = np.array(row[f"true_mhc_origin"])
    true_tcr_axes = np.array(row[f"true_tcr_axes"])
    true_tcr_origin = np.array(row[f"true_tcr_origin"])
    true_dgeom = np.array(row[f"true_dgeom"])

    true_u = align_mhc_to_origin(true_u, true_mhc_axes, true_mhc_origin)

    pred_cdr_coms = []
    true_cdr_coms = []

    for segid, chain_name, chain_species in zip(
        ["D", "E"], ["alpha", "beta"], ["tcr_1_species", "tcr_2_species"]
    ):

        tcr_pred_res = pred_u.select_atoms(f"segid {segid} and name CA").residues
        tcr_true_res = true_u.select_atoms(f"segid {segid} and name CA").residues

        tcr_true_res, tcr_pred_res = seq_align_residue_groups(
            tcr_true_res, tcr_pred_res
        )
        _, imgt_num, (start, stop) = annotate_tcr(
            tcr_true_res.sequence(format="string"),
            tcr_true_res.resindices,
            chain_name,
            # this is a best guess, anarci will select the correct organism
            # if this fails
            chain_species,
        )

        tcr_true_res = tcr_true_res[start:stop]
        tcr_pred_res = tcr_pred_res[start:stop]

        for cdr_region in [(27, 38), (56, 65), (81, 86), (104, 118)]:

            # 1,2,3: https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
            # 2.5: https://pmc.ncbi.nlm.nih.gov/articles/PMC5616171/
            tcr_imgt_num_subsel_other_cdr = list(
                np.argwhere(
                    ((imgt_num >= cdr_region[0]) & (imgt_num <= cdr_region[1]))
                )[:, 0]
            )

            true_cdr_coms.append(
                tcr_true_res[tcr_imgt_num_subsel_other_cdr]
                .atoms.select_atoms("name CA and (altloc A or not altloc [!?])")
                .center_of_mass()
            )

            pred_cdr_coms.append(
                tcr_pred_res[tcr_imgt_num_subsel_other_cdr]
                .atoms.select_atoms("name CA and (altloc A or not altloc [!?])")
                .center_of_mass()
            )

    true_cdr_COMs = np.vstack(true_cdr_coms)
    pred_cdr_COMs = np.vstack(pred_cdr_coms)

    weights = [1, 1, 1, 3, 1, 1, 1, 3]

    docking_rmsd = rms.rmsd(
        true_cdr_COMs,
        pred_cdr_COMs,
        weights=weights,
        center=False,
        superposition=False,
    )

    row[f"docking_rmsd_{inf_type}_{rank}"] = docking_rmsd
    return row


def true_pred_cdr_rmsd(row, inf_type, clean_pdb_dir, inf_dir, rank=None):
    sample_id = row[f"pred_sample_rank_{rank}"]

    if inf_type == "af3":
        inf_cls = AF3Output
        sample_kwarg = {"seed": sample_id}
    elif inf_type == "boltz":
        inf_cls = BoltzOutput
        sample_kwarg = {"sample_num": sample_id}
    else:
        raise ValueError

    pred_u = inf_cls(inf_dir / row["job_name"]).get_mda_universe(**sample_kwarg)
    pred_mhc_axes = np.array(row[f"pred_mhc_axes_{rank}"])
    pred_mhc_origin = np.array(row[f"pred_mhc_origin_{rank}"])
    pred_tcr_axes = np.array(row[f"pred_tcr_axes_{rank}"])
    pred_tcr_origin = np.array(row[f"pred_tcr_origin_{rank}"])
    pred_dgeom = np.array(row[f"pred_dgeom_{rank}"])

    pred_u = align_mhc_to_origin(pred_u, pred_mhc_axes, pred_mhc_origin)

    true_u = mda.Universe(clean_pdb_dir / f"{row['pdb']}.pdb")
    true_mhc_axes = np.array(row[f"true_mhc_axes"])
    true_mhc_origin = np.array(row[f"true_mhc_origin"])
    true_tcr_axes = np.array(row[f"true_tcr_axes"])
    true_tcr_origin = np.array(row[f"true_tcr_origin"])
    true_dgeom = np.array(row[f"true_dgeom"])

    true_u = align_mhc_to_origin(true_u, true_mhc_axes, true_mhc_origin)

    pred_positions = []
    true_positions = []
    weights = []

    for segid, chain_name, chain_species in zip(
        ["D", "E"], ["alpha", "beta"], ["tcr_1_species", "tcr_2_species"]
    ):

        tcr_pred_res = pred_u.select_atoms(f"segid {segid} and name CA").residues
        tcr_true_res = true_u.select_atoms(f"segid {segid} and name CA").residues

        tcr_true_res, tcr_pred_res = seq_align_residue_groups(
            tcr_true_res, tcr_pred_res
        )
        _, imgt_num, (start, stop) = annotate_tcr(
            tcr_true_res.sequence(format="string"),
            tcr_true_res.resindices,
            chain_name,
            # this is a best guess, anarci will select the correct organism
            # if this fails
            chain_species,
        )

        tcr_true_res = tcr_true_res[start:stop]
        tcr_pred_res = tcr_pred_res[start:stop]

        # 1,2,3: https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
        # 2.5: https://pmc.ncbi.nlm.nih.gov/articles/PMC5616171/
        tcr_imgt_num_subsel_other_cdr = list(
            np.argwhere(
                # 1
                ((imgt_num >= 27) & (imgt_num <= 38))
                # 2
                | ((imgt_num >= 56) & (imgt_num <= 65))
                # 2.5
                | ((imgt_num >= 81) & (imgt_num <= 86))
            )[:, 0]
        )

        # incl. conserved C and F
        tcr_imgt_num_subsel_cdr3 = list(
            np.argwhere((imgt_num >= 104) & (imgt_num <= 118))[:, 0]
        )

        weights_chain = [1 for i in range(len(tcr_imgt_num_subsel_other_cdr))] + [
            3 for i in range(len(tcr_imgt_num_subsel_cdr3))
        ]

        tcr_cdr_subsel = tcr_imgt_num_subsel_other_cdr + tcr_imgt_num_subsel_cdr3

        tcr_pdb_cdr_res = tcr_true_res[tcr_cdr_subsel]
        tcr_pred_cdr_res = tcr_pred_res[tcr_cdr_subsel]

        true_positions.append(
            tcr_pdb_cdr_res.atoms.select_atoms(
                "name CA and (altloc A or not altloc [!?])"
            ).positions
        )
        pred_positions.append(tcr_pred_cdr_res.atoms.select_atoms("name CA").positions)

        weights += weights_chain

    cdr_rmsd = rms.rmsd(
        np.concat(true_positions),
        np.concat(pred_positions),
        weights=weights,
        center=False,
        superposition=False,
    )

    row[f"cdr_rmsd_{inf_type}_{rank}"] = cdr_rmsd
    return row


def true_pred_peptide_rmsd(row, inf_type, clean_pdb_dir, inf_dir, rank=None):
    sample_id = row[f"pred_sample_rank_{rank}"]

    if inf_type == "af3":
        inf_cls = AF3Output
        sample_kwarg = {"seed": sample_id}
    elif inf_type == "boltz":
        inf_cls = BoltzOutput
        sample_kwarg = {"sample_num": sample_id}
    else:
        raise ValueError

    pred_u = inf_cls(inf_dir / row["job_name"]).get_mda_universe(**sample_kwarg)
    pred_mhc_axes = np.array(row[f"pred_mhc_axes_{rank}"])
    pred_mhc_origin = np.array(row[f"pred_mhc_origin_{rank}"])
    pred_tcr_axes = np.array(row[f"pred_tcr_axes_{rank}"])
    pred_tcr_origin = np.array(row[f"pred_tcr_origin_{rank}"])
    pred_dgeom = np.array(row[f"pred_dgeom_{rank}"])

    pred_u = align_mhc_to_origin(pred_u, pred_mhc_axes, pred_mhc_origin)

    true_u = mda.Universe(clean_pdb_dir / f"{row['pdb']}.pdb")
    true_mhc_axes = np.array(row[f"true_mhc_axes"])
    true_mhc_origin = np.array(row[f"true_mhc_origin"])
    true_tcr_axes = np.array(row[f"true_tcr_axes"])
    true_tcr_origin = np.array(row[f"true_tcr_origin"])
    true_dgeom = np.array(row[f"true_dgeom"])

    true_u = align_mhc_to_origin(true_u, true_mhc_axes, true_mhc_origin)

    segid_sel = "segid A"

    true_peptide_res = true_u.select_atoms(f"{segid_sel} and name CA").residues

    pred_peptide_res = pred_u.select_atoms(f"{segid_sel} and name CA").residues

    true_peptide_res, pred_peptide_res = seq_align_residue_groups(
        true_peptide_res, pred_peptide_res
    )

    # can't use backbone, 5tez
    # if (
    #     pred_peptide_res.atoms.select_atoms("backbone").positions.shape
    #     != true_peptide_res.atoms.select_atoms("backbone").positions.shape
    # ):
    #     raise ValueError(
    #         f"Peptide residue groups do not match in shape: {pred_peptide_res.shape} vs {true_peptide_res.shape}"
    #     )

    peptide_rmsd = rms.rmsd(
        pred_peptide_res.atoms.select_atoms("name CA").positions,
        true_peptide_res.atoms.select_atoms("name CA").positions,
        center=False,
        superposition=False,
    )

    row[f"peptide_rmsd_{inf_type}_{rank}"] = peptide_rmsd
    return row


def true_pred_mhc_rmsd(row, inf_type, clean_pdb_dir, inf_dir, rank=None):
    sample_id = row[f"pred_sample_rank_{rank}"]

    if inf_type == "af3":
        inf_cls = AF3Output
        sample_kwarg = {"seed": sample_id}
    elif inf_type == "boltz":
        inf_cls = BoltzOutput
        sample_kwarg = {"sample_num": sample_id}
    else:
        raise ValueError

    pred_u = inf_cls(inf_dir / row["job_name"]).get_mda_universe(**sample_kwarg)
    true_u = mda.Universe(clean_pdb_dir / f"{row['pdb']}.pdb")

    mhc_sels = ["(segid B) and name CA"]
    if row["mhc_class"] == "II":
        mhc_sels.append("(segid C) and name CA")

    pred_mhc_resgrp_ls = []
    true_mhc_resgrp_ls = []

    for mhc_sel in mhc_sels:
        pred_mhc_res = pred_u.select_atoms(mhc_sel).residues
        true_mhc_res = true_u.select_atoms(mhc_sel).residues

        true_mhc_res, pred_mhc_res = seq_align_residue_groups(
            true_mhc_res, pred_mhc_res
        )

        pred_mhc_resgrp_ls.append(pred_mhc_res)
        true_mhc_resgrp_ls.append(true_mhc_res)

    row[f"mhc_rmsd_{inf_type}_{rank}"] = rms.rmsd(
        sum(pred_mhc_resgrp_ls).atoms.select_atoms("backbone").positions,
        sum(true_mhc_resgrp_ls).atoms.select_atoms("backbone").positions,
        superposition=True,
    )
    return row


def true_pred_tcr_rmsd(row, inf_type, clean_pdb_dir, inf_dir, rank=None):
    sample_id = row[f"pred_sample_rank_{rank}"]

    if inf_type == "af3":
        inf_cls = AF3Output
        sample_kwarg = {"seed": sample_id}
    elif inf_type == "boltz":
        inf_cls = BoltzOutput
        sample_kwarg = {"sample_num": sample_id}
    else:
        raise ValueError

    pred_u = inf_cls(inf_dir / row["job_name"]).get_mda_universe(**sample_kwarg)
    true_u = mda.Universe(clean_pdb_dir / f"{row['pdb']}.pdb")

    tcr_sels = ["(segid D) and name CA", "(segid E) and name CA"]
    true_tcr_resgrp_ls = []
    pred_tcr_resgrp_ls = []

    for tcr_sel in tcr_sels:
        pred_tcr_res = pred_u.select_atoms(tcr_sel).residues
        true_tcr_res = true_u.select_atoms(tcr_sel).residues

        true_tcr_res, pred_tcr_res = seq_align_residue_groups(
            true_tcr_res, pred_tcr_res
        )

        pred_tcr_resgrp_ls.append(pred_tcr_res)
        true_tcr_resgrp_ls.append(true_tcr_res)

    row[f"tcr_rmsd_{inf_type}_{rank}"] = rms.rmsd(
        sum(true_tcr_resgrp_ls).atoms.select_atoms("backbone").positions,
        sum(pred_tcr_resgrp_ls).atoms.select_atoms("backbone").positions,
        superposition=True,
    )
    return row


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_parquet",
        type=str,
    )
    parser.add_argument("--cleaned_pdbs", type=str)
    parser.add_argument("--inference_type")
    parser.add_argument("--inference_dir", type=str)
    parser.add_argument("--true_tcrdock_pq", type=str)
    parser.add_argument("--pred_tcrdock_pq", type=str)
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
    )
    args = parser.parse_args()

    pdb_df = pl.read_parquet(args.input_parquet)

    true_tcrdock_pq = pl.read_parquet(args.true_tcrdock_pq)
    pred_tcrdock_pq = pl.read_parquet(args.pred_tcrdock_pq)

    pdb_df = pdb_df.join(
        true_tcrdock_pq.select(cs.starts_with("job_name") | cs.starts_with("true")),
        on="job_name",
    ).join(
        pred_tcrdock_pq.select(cs.starts_with("job_name") | cs.starts_with("pred")),
        on="job_name",
    )

    inf_dir = Path(args.inference_dir)
    cleaned_pdb_dir = Path(args.cleaned_pdbs)

    for rank in range(0, 5):
        pdb_df = split_apply_combine(
            pdb_df,
            true_pred_docking_rmsd,
            args.inference_type,
            cleaned_pdb_dir,
            inf_dir,
            rank=rank,
            chunksize=15,
        )

        pdb_df = split_apply_combine(
            pdb_df,
            true_pred_cdr_rmsd,
            args.inference_type,
            cleaned_pdb_dir,
            inf_dir,
            rank=rank,
            chunksize=15,
        )

        pdb_df = split_apply_combine(
            pdb_df,
            true_pred_peptide_rmsd,
            args.inference_type,
            cleaned_pdb_dir,
            inf_dir,
            rank=rank,
            chunksize=15,
        )

        pdb_df = split_apply_combine(
            pdb_df,
            true_pred_mhc_rmsd,
            args.inference_type,
            cleaned_pdb_dir,
            inf_dir,
            rank=rank,
            chunksize=15,
        )

        pdb_df = split_apply_combine(
            pdb_df,
            true_pred_tcr_rmsd,
            args.inference_type,
            cleaned_pdb_dir,
            inf_dir,
            rank=rank,
            chunksize=15,
        )

    pdb_df.write_parquet(args.output_path)
