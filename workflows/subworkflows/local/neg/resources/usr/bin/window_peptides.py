#!/usr/bin/env python
from tcrtrifold.utils import generate_job_name, FORMAT_ANTIGEN_COLS
import polars as pl

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--neg_df",
        type=str,
        help="Path to the negative triad dataframe",
    )
    parser.add_argument(
        "--output_pmhc_path",
        type=str,
    )
    parser.add_argument(
        "--output_triad_path",
        type=str,
    )
    parser.add_argument(
        "--window_size",
        type=int,
        default=9,
        help="Size of the peptide window to generate",
    )
    args = parser.parse_args()

    neg_df = pl.read_parquet(args.neg_df)

    antigen_with_id = generate_job_name(
        neg_df.select(FORMAT_ANTIGEN_COLS).unique(),
        [
            "peptide",
            "mhc_1_seq",
            "mhc_2_seq",
        ],
        name="antigen_name",
    ).with_columns(pl.col("peptide").alias("orig_peptide"))

    neg_df = neg_df.join(
        antigen_with_id,
        on=FORMAT_ANTIGEN_COLS,
    )

    neg_df = neg_df.with_columns(
        pl.col("job_name").alias("triad_name"),
    ).select(pl.exclude("job_name"))

    neg_df = neg_df.with_columns(
        pl.col("peptide")
        .map_elements(
            lambda x: [
                x[i : i + args.window_size]
                for i in range(0, len(x) - (args.window_size - 1))
            ],
            return_dtype=pl.List(str),
        )
        .alias("peptide")
    ).explode("peptide")

    neg_df = generate_job_name(
        neg_df,
        [
            "peptide",
            "mhc_1_seq",
            "mhc_2_seq",
            "tcr_1_seq",
            "tcr_2_seq",
        ],
    )

    win_antigen = neg_df.select(
        FORMAT_ANTIGEN_COLS + ["antigen_name", "orig_peptide"]
    ).unique()

    win_antigen = generate_job_name(
        win_antigen,
        [
            "peptide",
            "mhc_1_seq",
            "mhc_2_seq",
        ],
    )

    win_antigen.write_parquet(args.output_pmhc_path)

    neg_df.write_parquet(args.output_triad_path)
