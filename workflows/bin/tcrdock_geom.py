#!/usr/bin/env python

import tcrdock
from mdaf3.FeatureExtraction import serial_apply
from mdaf3.AF3OutputParser import AF3Output
import MDAnalysis as mda
from pathlib import Path
import polars as pl
import argparse
from collections import OrderedDict, defaultdict


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
            else:
                raise ValueError

            res_dict[fmt_atom_name] = pos

        pose_coords[res_key] = res_dict

    pose_ca_coords = ag.residues.atoms.select_atoms("name CA").positions
    pose_chains = ag.segments.segids.tolist()
    pose_chainseq = "/".join(
        [seg.residues.sequence(format="string") for seg in ag.segments]
    )
    pose_chainbounds = [int(seg.residues[0].ix) for seg in ag.segments] + [len(seq)]

    pose["resids"] = pose_resids
    pose["coords"] = pose_coords
    pose["sequence"] = seq
    pose["ca_coords"] = pose_ca_coords
    pose["chains"] = pose_chains
    pose["chainseq"] = pose_chainseq
    pose["chainbounds"] = pose_chainbounds

    return pose


def reorder_pose(pose, segid_map):
    """
    Reorder and rename segid chains in a pose dict.

    Args:
        pose (dict): original pose with:
            - "resids": list of (segid, residx) tuples
            - "coords": dict mapping (segid, residx) to {atom: np.array}
            - "sequence": concatenated amino acid string for all chains
        segid_map (dict): mapping old segid -> new segid

    Returns:
        dict: new pose with updated segids, renumbered residues, and reordered sequence
    """
    segid_to_entries = defaultdict(list)
    coords = pose["coords"]
    seq_iter = iter(pose["sequence"])

    for segid, residx in pose["resids"]:
        aa = next(seq_iter)
        atoms = coords[(segid, residx)]
        new_segid = segid_map.get(segid, segid)
        segid_to_entries[new_segid].append((residx, aa, atoms))

    new_resids = []
    new_coords = OrderedDict()
    new_sequence = ""

    for new_segid in segid_map.values():
        if new_segid not in segid_to_entries:
            continue
        entries = segid_to_entries[new_segid]
        for new_residx, (old_residx, aa, atoms) in enumerate(entries):
            key = (new_segid, str(new_residx))
            new_resids.append(key)
            new_coords[key] = atoms
            new_sequence += aa

    return {"resids": new_resids, "coords": new_coords, "sequence": new_sequence}


def extract_tcrdock_geom(row, af3_parent_dir, inf_key="job_name"):
    u = AF3Output(af3_parent_dir / row[inf_key]).get_mda_universe()

    pose = pose_from_mda_universe(u)
    fname = row[inf_key]
    mhc_class = 1 if row["mhc_class"] == "I" else "II"
    organism = row["mhc_1_species"]

    # the rest is modified from TCRdock/parse_tcr_pmhc_pdbfile.py
    num_chains = len(pose["chains"])
    if mhc_class == 1:
        if num_chains == 5:
            # remove B2M
            print(
                f"removing chain 2 from a 5-chain MHC class I pose; residue numbers "
                "in parsing output will not include this chain"
            )
            pose = tcrdock.pdblite.delete_chains(pose, [1])  # 0-indexed chain number
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

    dgeom = tcrdock.docking_geometry.DockingGeometry().from_stubs(mhc_stub, tcr_stub)

    return row


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_parquet_path",
        type=str,
    )
    parser.add_argument(
        "-d",
        "--inference_dir",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output_parquet_path",
        type=str,
    )
    args = parser.parse_args()

    df = pl.read_parquet(args.input_parquet_path)

    df = serial_apply(df, extract_tcrdock_geom, Path(args.inference_dir))

    df.write_parquet(
        args.output_parquet_path,
    )


from importlib import resources

example_pdb = (
    resources.files("tcrdock.db.pdb.ternary") / "1ao7.pdb.human.MH1.A-02.A.C.DE.pdb"
)
