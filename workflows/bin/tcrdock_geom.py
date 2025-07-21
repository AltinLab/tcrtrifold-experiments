#!/usr/bin/env python
import tcrdock
from mdaf3.FeatureExtraction import serial_apply, split_apply_combine
from mdaf3.AF3OutputParser import AF3Output
from mdaf3.BoltzOutputParser import BoltzOutput
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import polars as pl
import argparse
import re
from collections import OrderedDict


def get_best_seed_list_of_n_af3(af3_output):
    """
    Given the seeds of a compressed AF3 job, return a list "best_seeds"
    which gives the best seed given sequential runs of seeds [1..n]

    For example, given an AF3 job ran with 2 seeds (1, 2), where the first seed (1)
    performed the best, this method returns:

    [1, 1]

    The purpose of this is to simulate the effect of having run a smaller number of seeds
    on downstream metrics
    """

    with af3_output._get_h5_handle() as hf:
        ranking_score_np = hf["ranking_scores"]["ranking_score"][:]
        seed_np = hf["ranking_scores"]["seed"][:]

    best_seeds = []

    for seed in seed_np:
        focal_seed_idx = np.argwhere(seed_np <= seed)
        best_focal_seed_idx = np.argmax(ranking_score_np[focal_seed_idx])
        best_seeds.append(int(seed_np[focal_seed_idx[best_focal_seed_idx]][0]))

    return best_seeds


# def get_best_seed_list_of_n_boltz(boltz_output):
#     """
#     Performs the same function as `get_best_seed_list_of_n_af3` but with the boltz
#     output format. Boltz does not have seeds, but has diffusion samples with greater temperature.
#     """
#     dir = boltz_output.dir_path
#     for sample_conf in dir.glob("confidence_*.json"):
#         seed_np


def get_ascending_sample_list_boltz(boltz_output):
    """
    Performs the same function as `get_best_seed_list_of_n_af3` but with the boltz
    output format. Boltz does not have seeds, but has diffusion samples with greater temperature.
    """
    dir = boltz_output.dir_path
    sample_nums = []
    for sample_conf in dir.glob("confidence_*.json"):
        sample_num = re.findall(r"\d+", sample_conf.stem)[-1]
        sample_nums.append(sample_num)

    return sorted(sample_nums, reverse=True)


def reorder_universe(universe, reorder_map):

    chain_us = []

    reorder_map_srt = OrderedDict(sorted(reorder_map.items(), key=lambda x: x[1]))

    for from_segid, to_segid in reorder_map_srt.items():

        if from_segid not in universe.segments.segids:
            raise ValueError

        if to_segid not in universe.segments.segids:
            raise ValueError

        chain_u = mda.Merge(universe.select_atoms(f"segid {from_segid}"))
        chain_u.segments.segids = to_segid
        chain_u.atoms.chainIDs = [to_segid] * len(chain_u.atoms)
        chain_us.append(chain_u.atoms)

    return mda.Merge(*chain_us)


def pose_from_mda_universe(u):
    """
    A pose is a dict object for use with the TCRdock package:

    {
        "resids" : [
            (<segid>, str(<residx>))
        ],
        "coords : {
            (<segid>, str(<residx>)) : {
                <atom name> : np.array(x,y,z),
            },
        },
        "sequence" : <AA seq>
        "ca_coords": np.array[
            [x,y,z],
            ...
        ],
        "chains": [
            <segid>,
            ...
        ],
        "chainseq": "<seg 1 seq>/<seg 2 seq>/...",
        "chainbounds": [
            <offset of seg 1 seq in sequence>,
            ...
            <length of sequence>,
        ],
    }

    In TCRdock, the order of chains is important. For class II, it is always:

    [mhc_1, mhc_2, peptide, tcr_1, tcr_2]

    For class I:

    [mhc_1, Optional[mhc_2], peptide, tcr_1, tcr_2]
    """
    pose = {}
    ag = u

    seq = ag.residues.sequence(format="string")
    segids = ag.residues.segids.tolist()
    residx = ag.residues.resindices.astype("str").tolist()

    pose_resids = list(zip(segids, residx))

    pose_coords = {}

    for res_key, res in zip(pose_resids, ag.residues):

        res_dict = OrderedDict()
        for atom_name, pos in zip(res.atoms.names.tolist(), res.atoms.positions):
            if len(atom_name) == 3:
                fmt_atom_name = atom_name + " "
            elif len(atom_name) == 2:
                fmt_atom_name = " " + atom_name + " "
            elif len(atom_name) == 1:
                fmt_atom_name = " " + atom_name + "  "
            elif len(atom_name) == 4:
                fmt_atom_name = atom_name
            else:
                raise ValueError

            res_dict[fmt_atom_name] = pos

        pose_coords[res_key] = res_dict

    # if we needed to extract derived data manually:
    # pose_ca_coords = ag.residues.atoms.select_atoms("name CA").positions
    # pose_chains = ag.segments.segids.tolist()
    # pose_chainseq = "/".join(
    #     [seg.residues.sequence(format="string") for seg in ag.segments]
    # )
    # pose_chainbounds = [int(seg.residues[0].ix) for seg in ag.segments] + [len(seq)]

    pose["resids"] = pose_resids
    pose["coords"] = pose_coords
    pose["sequence"] = seq

    pose = tcrdock.pdblite.update_derived_data(pose)
    return pose


def extract_tcrdock_geom(
    row, topology_parent_dir, from_true_struct=False, inf_type=None
):

    if from_true_struct:
        # u = mda.Universe(topology_parent_dir / f"{row['pdb']}.pdb")
        fname = row["pdb"]
        samples = [None]
    else:
        if inf_type == "af3":
            inf_output = AF3Output(topology_parent_dir / row["job_name"])
            samples = get_best_seed_list_of_n_af3(inf_output)
        else:
            inf_output = BoltzOutput(topology_parent_dir / row["job_name"])
            # can't actually detemine "best seed of n" for boltz
            # so just give ascending list of samples based on quality
            samples = get_ascending_sample_list_boltz(inf_output)

        fname = row["job_name"]

    for i, sample in enumerate(samples):

        if from_true_struct:
            u = mda.Universe(topology_parent_dir / f"{row['pdb']}.pdb")

        else:

            if inf_type == "af3":
                u = inf_output.get_mda_universe(seed=sample)

            elif inf_type == "boltz":
                u = inf_output.get_mda_universe(sample_num=sample)

        if row["mhc_class"] == "I" and row["mhc_2_seq"] is None:
            u = reorder_universe(u, {"A": "B", "B": "A", "D": "D", "E": "E"})

        else:
            u = reorder_universe(u, {"A": "C", "B": "A", "C": "B", "D": "D", "E": "E"})

        pose = pose_from_mda_universe(u)
        mhc_class = 1 if row["mhc_class"] == "I" else 2
        organism = row["mhc_1_species"]

        # the rest is modified from TCRdock/parse_tcr_pmhc_pdbfile.py
        num_chains = len(pose["chains"])
        if mhc_class == 1:
            if num_chains == 5:
                # remove B2M
                print(
                    "removing chain 2 from a 5-chain MHC class I pose; residue numbers "
                    "in parsing output will not include this chain"
                )
                pose = tcrdock.pdblite.delete_chains(
                    pose, [1]
                )  # 0-indexed chain number
                num_chains = len(pose["chains"])
            else:
                assert (
                    num_chains == 4
                ), f"MHC-I pdbfile {fname} should have 4 or 5 chains, see --help message"
            cs = pose["chainseq"].split("/")
            mhc_aseq, pep_seq, tcr_aseq, tcr_bseq = cs
            mhc_bseq = None
        else:
            assert (
                num_chains == 5
            ), f"MHC-II pdbfile {fname} should have 5 chains, see --help message"
            cs = pose["chainseq"].split("/")
            mhc_aseq, mhc_bseq, pep_seq, tcr_aseq, tcr_bseq = cs

        tdinfo = tcrdock.tcrdock_info.TCRdockInfo().from_sequences(
            organism,
            mhc_class,
            mhc_aseq,
            mhc_bseq,
            pep_seq,
            tcr_aseq,
            tcr_bseq,
        )

        # these are the MHC and TCR reference frames (aka 'stubs')
        mhc_stub = tcrdock.mhc_util.get_mhc_stub(pose, tdinfo)
        tcr_stub = tcrdock.tcr_util.get_tcr_stub(pose, tdinfo)

        dgeom = tcrdock.docking_geometry.DockingGeometry().from_stubs(
            mhc_stub, tcr_stub
        )

        ## further modifications from original script ##
        # store everything in dicts
        # conv np.ndarray to list so it plays nice with polars
        mhc_stub_dict = {k: v.tolist() for k, v in mhc_stub.items()}
        tcr_stub_dict = {k: v.tolist() for k, v in tcr_stub.items()}

        dgeom_dict = {
            "d": float(dgeom.d),
            "mhc_unit_x_is_negative": bool(dgeom.mhc_unit_x_is_negative),
            "mhc_unit_y": float(dgeom.mhc_unit_y),
            "mhc_unit_z": float(dgeom.mhc_unit_z),
            "tcr_unit_x_is_negative": bool(dgeom.tcr_unit_x_is_negative),
            "tcr_unit_y": float(dgeom.tcr_unit_y),
            "tcr_unit_z": float(dgeom.tcr_unit_z),
            "torsion": float(dgeom.torsion),
        }

        if from_true_struct:
            row["true_mhc_axes"] = mhc_stub_dict["axes"]
            row["true_mhc_origin"] = mhc_stub_dict["origin"]
            row["true_tcr_axes"] = tcr_stub_dict["axes"]
            row["true_tcr_origin"] = tcr_stub_dict["origin"]
            row["true_dgeom"] = dgeom_dict

        else:
            row[f"pred_mhc_axes_{i}"] = mhc_stub_dict["axes"]
            row[f"pred_mhc_origin_{i}"] = mhc_stub_dict["origin"]
            row[f"pred_tcr_axes_{i}"] = tcr_stub_dict["axes"]
            row[f"pred_tcr_origin_{i}"] = tcr_stub_dict["origin"]
            row[f"pred_dgeom_{i}"] = dgeom_dict
            row[f"pred_sample_rank_{i}"] = sample

    return row


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_parquet",
        type=str,
    )
    parser.add_argument(
        "--topology_path",
        type=str,
    )
    parser.add_argument("--from_true_struct", action="store_true")
    parser.add_argument("--inference_type")
    parser.add_argument(
        "-o",
        "--output_parquet_path",
        type=str,
    )
    args = parser.parse_args()

    if not args.from_true_struct:
        inf_type = args.inference_type
    else:
        inf_type = None

    df = pl.read_parquet(args.input_parquet)

    df = split_apply_combine(
        df,
        extract_tcrdock_geom,
        Path(args.topology_path),
        args.from_true_struct,
        inf_type=inf_type,
        chunksize=15,
    )

    df.write_parquet(
        args.output_parquet_path,
    )
