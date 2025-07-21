#!/usr/bin/env python
"""
Generate negative triads given a base dataframe and a list of "supp_dfs" dataframes
This process is deterministic- polars' random seed is set, and all
unique/group_by/join/partition_by operations are set to maintain order.
"""
from tcrtrifold.neg_creation import (
    generate_all_possible_negs,
    sample_to,
)
from tcrtrifold.utils import FORMAT_ANTIGEN_COLS
import argparse
import polars as pl
import random

random.seed(42)
pl.set_random_seed(42)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base_df",
        type=str,
    )
    parser.add_argument("--supp_dfs", type=str)
    parser.add_argument("--neg_depth", type=int)
    parser.add_argument(
        "-d",
        "--discard_path",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
    )
    args = parser.parse_args()

    base_df = pl.read_parquet(args.base_df).sort(by="job_name")

    # for the sake of this script, set null mhc_2_seqs to ""
    base_df = base_df.with_columns(
        pl.when(pl.col("mhc_2_seq").is_null())
        .then(pl.lit(""))
        .otherwise(pl.col("mhc_2_seq"))
        .alias("mhc_2_seq"),
        pl.when(pl.col("mhc_2_name").is_null())
        .then(pl.lit(""))
        .otherwise(pl.col("mhc_2_name"))
        .alias("mhc_2_name"),
        pl.when(pl.col("mhc_2_species").is_null())
        .then(pl.lit(""))
        .otherwise(pl.col("mhc_2_species"))
        .alias("mhc_2_species"),
        pl.when(pl.col("mhc_2_chain").is_null())
        .then(pl.lit(""))
        .otherwise(pl.col("mhc_2_chain"))
        .alias("mhc_2_chain"),
    )

    antigen_df = base_df.group_by(FORMAT_ANTIGEN_COLS, maintain_order=True).agg(
        (pl.len() * args.neg_depth).alias("needed_negs")
    )

    supp_dfs = [
        pl.read_parquet(supp_df_path).sort(by="job_name")
        for supp_df_path in args.supp_dfs.split(",")
    ]

    remaining_antigens = antigen_df
    neg_dfs = []

    for supp_df in supp_dfs:
        all_negs_df = generate_all_possible_negs(base_df, supp_df)
        negs_df, remaining_antigens = sample_to(remaining_antigens, all_negs_df)

        neg_dfs.append(negs_df)

        if remaining_antigens is None:
            break

    out_df = pl.concat([base_df] + neg_dfs, how="diagonal_relaxed")

    out_df = out_df.with_columns(
        pl.when(pl.col("mhc_2_seq") == "")
        .then(pl.lit(None))
        .otherwise(pl.col("mhc_2_seq"))
        .alias("mhc_2_seq"),
        pl.when(pl.col("mhc_2_name") == "")
        .then(pl.lit(None))
        .otherwise(pl.col("mhc_2_name"))
        .alias("mhc_2_name"),
        pl.when(pl.col("mhc_2_species") == "")
        .then(pl.lit(None))
        .otherwise(pl.col("mhc_2_species"))
        .alias("mhc_2_species"),
        pl.when(pl.col("mhc_2_chain") == "")
        .then(pl.lit(None))
        .otherwise(pl.col("mhc_2_chain"))
        .alias("mhc_2_chain"),
    )

    out_df.write_parquet(args.output_path)

    if remaining_antigens is not None:
        remaining_antigens.write_parquet(args.discard_path)
