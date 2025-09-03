#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
import hashlib


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


def get_msa(msa_cache_dir, protein_type, seq):

    h = hashlib.sha256()

    if protein_type == "peptide":
        fname = seq + ".json"
    else:
        h.update(seq.encode("utf-8"))
        fname = h.hexdigest() + ".json"

    msa_path = msa_cache_dir / protein_type / fname

    if msa_path.is_file():
        with open(msa_path) as f:
            msa = json.load(f)
    else:
        raise FileNotFoundError(
            f"MSA file not found: {msa_path}. Please ensure the MSA cache directory is correct."
        )

    return msa


def main():
    parser = argparse.ArgumentParser(
        description="Compose Alphafold3 input JSON by querying VAST for chain MSA information."
    )
    parser.add_argument("--job_name", type=str, required=True, help="Job name")
    parser.add_argument(
        "--fasta",
        type=str,
        required=True,
        help="Path to fasta file",
    )
    parser.add_argument(
        "--skip_msa",
        type=str,
        help="Skip MSA for sequence at index i",
    )
    parser.add_argument(
        "--protein_types",
        type=str,
        required=True,
        help="Comma separated list of protein types",
    )
    parser.add_argument(
        "--msa_cache_dir",
        type=str,
        required=True,
        help="Directory to retrieve MSAs from",
    )
    parser.add_argument(
        "--inf_dir",
        type=str,
        required=False,
        help="Directory to check for existing inference results",
    )
    parser.add_argument(
        "--segids",
        type=str,
        required=False,
        help="Comma separated list of segids (chain IDs) the same length as the number of proteins",
    )
    parser.add_argument(
        "--seeds",
        type=str,
        required=False,
        default="42",
        help="Comma separated list of model seeds",
    )
    parser.add_argument(
        "--check_inf_exists",
        action="store_true",
        help="Check if inference already exists in the specified directory",
    )

    args = parser.parse_args()

    if args.skip_msa:
        skip_msa = set([int(i) for i in args.skip_msa.split(",")])

    print(f"Checking inference directory: {args.inf_dir}")

    if args.check_inf_exists:
        inf_dir = Path(args.inf_dir)

        if (inf_dir / args.job_name).is_dir():

            print(
                f"Skipping job {args.job_name} as inference already exists in {inf_dir}."
            )
            return

        else:
            print(f"No existing inference found for job {args.job_name} in {inf_dir}.")

    segids = args.segids.split(",") if args.segids else None

    protein_type = list(args.protein_types.split(","))

    seeds = [int(seed) for seed in args.seeds.split(",")]

    msa_cache_dir = Path(args.msa_cache_dir)

    msas = []

    seqs = read_fasta_seqs(args.fasta)

    if segids is not None and len(segids) != len(seqs):
        raise ValueError

    for i, seq in enumerate(seqs):

        if segids is not None:
            segid = segids[i]
        else:
            segid = chr(65 + 1)

        if args.skip_msa and i in skip_msa:
            msa = {
                "protein": {
                    "id": segid,
                    "sequence": seq,
                    "unpairedMsa": f">query\n{seq}\n",
                    "pairedMsa": f">query\n{seq}\n",
                    "templates": [],
                },
            }
        else:
            msa = get_msa(msa_cache_dir, protein_type[i], seq)
            msa["id"] = segid
            msa = {"protein": msa}
        msas.append(msa)

    final_json = {
        "name": args.job_name,
        "modelSeeds": seeds,
        "sequences": msas,
        "dialect": "alphafold3",
        "version": 1,
    }
    with open(args.job_name + ".json", "w") as f:
        json.dump(final_json, f, indent=2)


if __name__ == "__main__":
    main()
