#!/usr/bin/env python
from tcrtrifold.utils import (
    generate_job_name,
    FORMAT_COLS,
    FORMAT_ANTIGEN_COLS,
    TCRDIST_COLS,
)
from tcrtrifold.tcr import (
    extract_tcrdist_cols,
)
from mdaf3.FeatureExtraction import serial_apply
import polars as pl
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--raw_csv_path",
        type=str,
    )
    parser.add_argument(
        "-op",
        "--output_pmhc_path",
        type=str,
    )
    parser.add_argument(
        "-ot",
        "--output_triad_path",
        type=str,
    )
    args = parser.parse_args()

    addtl_cols = TCRDIST_COLS

    cresta = pl.read_csv(args.raw_csv_path)

    cresta = serial_apply(
        cresta,
        extract_tcrdist_cols,
    )

    cresta.select(FORMAT_COLS + addtl_cols).write_parquet(
        args.output_triad_path,
    )

    cresta_antigen = cresta.select(FORMAT_ANTIGEN_COLS).unique()
    cresta_antigen = generate_job_name(
        cresta_antigen,
        ["peptide", "mhc_1_seq", "mhc_2_seq"],
    )

    cresta_antigen.select(["job_name"] + FORMAT_ANTIGEN_COLS).write_parquet(
        args.output_pmhc_path,
    )
