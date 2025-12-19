#!/usr/bin/env python3

import argparse
import yaml
from pathlib import Path
import hashlib
import json


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
    parser = argparse.ArgumentParser(description="Compose Boltz input YAML.")
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
    parser.add_argument(
        "--inf_dir",
        type=str,
        required=False,
        help="Directory where inference files will be stored",
    )
    parser.add_argument(
        "--msa_cache_dir",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    seqs = read_fasta_seqs(args.fasta_path)

    if args.skip_msa:
        skip_msa = set([int(i) for i in args.skip_msa.split(",")])

    inf_dir = Path(args.inf_dir)
    msa_cache_dir = Path(args.msa_cache_dir)

    if (inf_dir / args.job_name).is_dir():

        print(f"Skipping job {args.job_name} as inference already exists in {inf_dir}.")
        return

    segids = args.segids.split(",") if args.segids else None

    if len(seqs) == 4:
        protein_types = ["peptide", "mhc", "tcr", "tcr"]
    else:
        protein_types = ["peptide", "mhc", "mhc", "tcr", "tcr"]

    yaml_data = {"sequences": [], "templates": []}

    for i, seq in enumerate(seqs):

        if segids is not None:
            segid = segids[i]
        else:
            segid = chr(65 + 1)
        entity = {
            "protein": {
                "id": segid,
                "sequence": seq,
            }
        }
        if i in skip_msa:
            entity["protein"]["msa"] = "empty"

        yaml_data["sequences"].append(entity)

        if i != 0:
            af3_msa = get_msa(msa_cache_dir, protein_types[i], seq)
            chain_templates = af3_msa["templates"]

            for j, chain_template in enumerate(chain_templates):
                fname = f"{segid}_{j}.cif"
                yaml_data["templates"].append({"cif": fname, "chain_id": segid})
                cif_dat = chain_template["mmcif"]

                with open(fname, "w") as f:
                    f.write(cif_dat)

    with open(args.job_name + ".yaml", "w") as f:
        yaml.safe_dump(yaml_data, f)


if __name__ == "__main__":
    main()
