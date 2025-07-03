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

    out_df.write_parquet(args.output_path)

    if remaining_antigens is not None:
        remaining_antigens.write_parquet(args.discard_path)
