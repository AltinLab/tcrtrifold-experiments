import gzip
from alphafold3.structure.parsing import from_mmcif
from alphafold3.model.scoring.chirality import compare_chirality
import polars as pl
from pathlib import Path
import argparse


def read_gz_to_str(path: str, encoding: str = "utf-8") -> str:
    with gzip.open(path, mode="rt", encoding=encoding) as f:
        return f.read()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_parquet",
        type=str,
    )
    parser.add_argument(
        "--inf_dir",
        type=str,
    )
    parser.add_argument(
        "--output_dir",
        type=str,
    )
    args = parser.parse_args()

    inf_dir = Path(args.inf_dir)

    pq = pl.read_parquet(args.input_parquet)

    for row in pq.iter_rows(named=True):

        cif_str = read_gz_to_str(
            (inf_dir / row["job_name"]) / (row["job_name"] + "_model.cif.gz")
        )
        struct = from_mmcif(cif_str)

        chirality_dict = compare_chirality(struct)

        if chirality_dict:
            print("hello world")
