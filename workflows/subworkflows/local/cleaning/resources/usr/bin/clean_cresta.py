#!/usr/bin/env python
from tcrtrifold.utils import (
    generate_job_name,
    FORMAT_COLS,
    FORMAT_ANTIGEN_COLS,
    TCRDIST_COLS,
)
from tcrtrifold.mhc import (
    B2M_HUMAN_SEQ,
    HLACodeWebConverter,
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
    # parser.add_argument(
    #     "--peptide_correction_csv",
    #     type=str,
    # )
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

    addtl_cols = TCRDIST_COLS + ["suspected_9mer"]

    cresta = pl.read_csv(args.raw_csv_path)

    human_conv = HLACodeWebConverter()

    cresta = cresta.with_columns(
        pl.col("mhc_1_name")
        .map_elements(
            lambda x: human_conv.get_sequence(x, top_only=True),
            return_dtype=pl.String,
        )
        .alias("mhc_1_seq"),
        pl.col("mhc_2_name")
        .map_elements(
            lambda x: human_conv.get_sequence(x, top_only=True),
            return_dtype=pl.String,
        )
        .alias("mhc_2_seq"),
    )

    # peptides = pl.read_csv(args.peptide_correction_csv)

    # cresta = (
    #     cresta.join(peptides, left_on="peptide", right_on="orig_peptide")
    #     .select(pl.exclude("peptide"))
    #     .rename({"peptide_right": "peptide"})
    # )

    # overwrite job name, but rename it "triad_name"
    cresta = generate_job_name(
        cresta,
        [
            "peptide",
            "mhc_1_seq",
            "mhc_2_seq",
            "tcr_1_seq",
            "tcr_2_seq",
        ],
    )

    cresta = serial_apply(
        cresta,
        extract_tcrdist_cols,
    )

    cresta.select(FORMAT_COLS + addtl_cols).unique().write_parquet(
        args.output_triad_path,
    )

    cresta_antigen = cresta.select(FORMAT_ANTIGEN_COLS + ["suspected_9mer"]).unique()
    cresta_antigen = generate_job_name(
        cresta_antigen,
        ["peptide", "mhc_1_seq", "mhc_2_seq"],
    )

    cresta_antigen.select(
        ["job_name"] + FORMAT_ANTIGEN_COLS + ["suspected_9mer"]
    ).write_parquet(
        args.output_pmhc_path,
    )
