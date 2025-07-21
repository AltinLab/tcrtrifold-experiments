#!/usr/bin/env python3

import argparse
import yaml

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compose Boltz input YAML."
    )
    parser.add_argument("-jn", "--job_name", type=str, required=True, help="Job name")
    parser.add_argument(
        "-f",
        "--fasta_path",
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
        "--segids",
        type=str,
        required=False,
        help="Comma separated list of segids (chain IDs) the same length as the number of proteins",
    )

    args = parser.parse_args()

    seqs = read_fasta_seqs(args.fasta_path)

    if args.skip_msa:
        skip_msa = set([int(i) for i in args.skip_msa.split(",")])

    segids = args.segids.split(",") if args.segids else None

    yaml_data = {"sequences" : []}

    for i, seq in enumerate(seqs):

        if segids is not None:
            segid = segids[i]
        else:
            segid = chr(65 + 1)
        entity = {
            "protein": {
                "id": segid,
                "sequence" : seq,
            }
        }
        if i in skip_msa:
            entity["protein"]["msa"] = "empty"

        yaml_data["sequences"].append(entity)

    with open(args.job_name + ".yaml", "w") as f:
        yaml.safe_dump(yaml_data, f, sort_keys=False)