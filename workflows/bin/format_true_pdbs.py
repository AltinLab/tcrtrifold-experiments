#!/usr/bin/env python
import polars as pl
import requests
import MDAnalysis as mda
from fuzzysearch import find_near_matches
from pathlib import Path
from io import StringIO
from mdaf3.FeatureExtraction import serial_apply
import argparse
import warnings
import tempfile

_tmp_dir = tempfile.TemporaryDirectory(prefix="tcrtrifold_")


def download_pdb(pdb_id):

    r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb")
    fmt = "PDB"
    file_like = StringIO(r.text)
    try:
        r.raise_for_status()
    except Exception as e:
        r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.cif")
        fmt = "CIF"
        r.raise_for_status()

        fname = Path(_tmp_dir.name) / f"{pdb_id}.cif"
        with open(fname, "w") as f:
            f.write(r.text)
        file_like = fname.as_posix()

    return file_like, fmt


def align_residue_group_to_seq(rg, seq):
    """
    Given the curated fasta seq of a PDB ID, sample the
    MDAnalysis ResidueGroup that represents the structure
    down to the portion that aligns with the curated seq

    We do this because sometimes multiple chains were stored
    under one segment ID in PDB
    """
    # if the residue group is about the same length as the seq,
    # this likely isn't a multiple-chains-under-one-segid
    # situation
    if (len(seq) / len(rg)) >= 0.8:
        return rg

    # otherwise, the residue group contains more than just the seq
    # most common case is one chain in the RG containing peptide + MHC I
    m = find_near_matches(seq, rg.sequence(format="string"), max_l_dist=8)
    if not m:
        raise ValueError("ResidueGroup not found in seq")

    if len(m) > 1:
        # find best match
        m.sort(key=lambda x: x.dist)
        warnings.warn(f"Found multiple matches, choosing {m[0]}")

    match = m[0]

    return rg[match.start : match.end]


def tcrtrifold_fmt_pdb(row, output_path):
    in_mem_pdb, fmt = download_pdb(row["pdb"])

    u = mda.Universe(in_mem_pdb, topology_format=fmt)

    sub_univ = []

    for segid_colname, curated_fasta_seq_colname, rename_segid in zip(
        [
            "peptide_segid",
            "mhc_1_segid",
            "mhc_2_segid",
            "tcr_1_segid",
            "tcr_2_segid",
        ],
        ["peptide", "mhc_1_seq", "mhc_2_seq", "tcr_1_seq", "tcr_2_seq"],
        ["A", "B", "C", "D", "E"],
    ):
        segid = row[segid_colname]

        if segid is None:
            if segid_colname != "mhc_2_segid":
                raise ValueError
            else:
                continue

        subset_rg = align_residue_group_to_seq(
            u.select_atoms(f"segid {segid} and name CA and record_type ATOM").residues,
            row[curated_fasta_seq_colname],
        )

        chain_u = mda.Merge(
            subset_rg.atoms.select_atoms(
                "record_type ATOM and (altloc A or not altloc [!?])"
            )
        )
        chain_u.segments.segids = rename_segid
        chain_u.atoms.chainIDs = [rename_segid] * len(chain_u.atoms)
        sub_univ.append(chain_u.atoms)

    fmt_u = mda.Merge(*sub_univ)

    with mda.Writer(output_path / (row["pdb"] + ".pdb")) as W:

        W.write(fmt_u.atoms)

    # noop
    return row


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pdb_parquet",
        type=str,
    )
    parser.add_argument(
        "--output_dir",
        type=str,
    )
    args = parser.parse_args()

    pdb_df = pl.read_parquet(args.pdb_parquet)

    serial_apply(pdb_df, tcrtrifold_fmt_pdb, Path(args.output_dir))
