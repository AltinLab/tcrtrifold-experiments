#!/usr/bin/env python
from tcrtrifold.utils import generate_job_name
from tcrtrifold.tcr import shorten_tcr_to_vregion
import polars as pl
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--triad_parquet",
        type=str,
    )
    parser.add_argument(
        "--output_triad",
        type=str,
    )
    parser.add_argument(
        "--output_peptide",
        type=str,
    )
    parser.add_argument(
        "--output_mhc_1",
        type=str,
    )
    parser.add_argument(
        "--output_mhc_2",
        type=str,
    )
    parser.add_argument(
        "--output_tcr_1",
        type=str,
    )
    parser.add_argument(
        "--output_tcr_2",
        type=str,
    )
    args = parser.parse_args()

    triad = pl.read_parquet(args.triad_parquet)

    triad = generate_job_name(
        triad,
        ["peptide"],
        name="peptide_hash",
    )
    triad = generate_job_name(
        triad,
        ["mhc_1_seq"],
        name="mhc_1_hash",
    )
    triad = generate_job_name(
        triad,
        ["mhc_2_seq"],
        name="mhc_2_hash",
    )
    triad = generate_job_name(
        triad,
        ["tcr_1_seq"],
        name="tcr_1_hash",
    )
    triad = generate_job_name(
        triad,
        ["tcr_2_seq"],
        name="tcr_2_hash",
    )

    peptide = (
        triad.select(["peptide", "peptide_hash"])
        .unique()
        .rename({"peptide": "seq", "peptide_hash": "name"})
    )
    mhc_1 = (
        triad.select(
            [
                "mhc_1_seq",
                "mhc_1_hash",
            ]
        )
        .unique()
        .rename(
            {
                "mhc_1_seq": "seq",
                "mhc_1_hash": "name",
            }
        )
    )

    mhc_2 = (
        triad.select(
            [
                "mhc_2_seq",
                "mhc_2_hash",
            ]
        )
        .unique()
        .rename(
            {
                "mhc_2_seq": "seq",
                "mhc_2_hash": "name",
            }
        )
    )
    tcr_1 = (
        triad.select(
            [
                "tcr_1_seq",
                "tcr_1_hash",
            ]
        )
        .unique()
        .with_columns(
            pl.col("tcr_1_seq")
            .map_elements(
                lambda x: shorten_tcr_to_vregion(x, "alpha", "human"),
                return_dtype=pl.String,
            )
            .alias("tcr_1_seq")
        )
        .rename(
            {
                "tcr_1_seq": "seq",
                "tcr_1_hash": "name",
            }
        )
    )
    tcr_2 = (
        triad.select(
            [
                "tcr_2_seq",
                "tcr_2_hash",
            ]
        )
        .with_columns(
            pl.col("tcr_2_seq")
            .map_elements(
                lambda x: shorten_tcr_to_vregion(x, "beta", "human"),
                return_dtype=pl.String,
            )
            .alias("tcr_2_seq")
        )
        .unique()
        .rename(
            {
                "tcr_2_seq": "seq",
                "tcr_2_hash": "name",
            }
        )
    )

    triad.write_parquet(args.output_triad)
    peptide.write_parquet(args.output_peptide)
    mhc_1.write_parquet(args.output_mhc_1)
    mhc_2.write_parquet(args.output_mhc_2)
    tcr_1.write_parquet(args.output_tcr_1)
    tcr_2.write_parquet(args.output_tcr_2)
