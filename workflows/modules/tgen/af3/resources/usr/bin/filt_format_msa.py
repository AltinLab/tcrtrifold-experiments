#!/usr/bin/env python3
import argparse
import os
import json
from pathlib import Path


def read_fasta_seqs(path):
    """
    Read a FASTA file and return a list of sequences (strings),
    concatenating multi-line records correctly.
    """
    seqs = []
    current_seq = []

    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                # If we were in the middle of a sequence, save it.
                if current_seq:
                    seqs.append("".join(current_seq))
                    current_seq = []
                # (We skip the header itself; if you need headers, collect them here.)
            else:
                # Append this line to the current sequence buffer
                current_seq.append(line)

        # After the loop, make sure to save the last sequence
        if current_seq:
            seqs.append("".join(current_seq))

    return seqs


def is_msa_stored(msa_cache_dir, protein_type, job_name):

    msa_path = msa_cache_dir / protein_type / (job_name + ".json")

    if msa_path.is_file():
        return True
    else:
        print(
            f"MSA file not found: {msa_path}. Please ensure the MSA cache directory is correct."
        )
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check if an MSA entry is missing from SQLite DB."
    )
    parser.add_argument("--protein_type", type=str, required=True, help="Protein type")
    parser.add_argument("--fasta", type=str, required=True, help="Protein sequence")
    parser.add_argument(
        "--force", action="store_true", required=False, help="Force update MSA"
    )
    parser.add_argument(
        "--job_name",
        type=str,
        required=True,
        help="Job name for the MSA entry",
    )
    parser.add_argument(
        "--msa_cache_dir",
        type=str,
        required=True,
        help="Directory to retrieve MSAs from",
    )
    args = parser.parse_args()
    seq = read_fasta_seqs(args.fasta)[0]

    msa_cache_dir = Path(args.msa_cache_dir)

    if args.force or not is_msa_stored(msa_cache_dir, args.protein_type, args.job_name):

        with open(args.job_name + ".json", "w") as f:
            json_dict = {
                "name": "af3-single-chain-msa",
                "modelSeeds": [42],
                "sequences": [{"protein": {"id": "A", "sequence": seq}}],
                "dialect": "alphafold3",
                "version": 1,
            }
            json.dump(json_dict, f, indent=2)
