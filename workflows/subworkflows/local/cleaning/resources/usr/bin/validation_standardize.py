#!/usr/bin/env python
from tcrtrifold.tcr import shorten_both_tcrs
from tcrtrifold.utils import FORMAT_ANTIGEN_COLS, generate_job_name
from mdaf3.FeatureExtraction import serial_apply
import argparse
import polars as pl


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pdb_parquet",
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

    val = (
        pl.read_parquet(args.pdb_parquet)
        .filter(~pl.col("replication"), pl.col("mhc_class") == "I")
        .select(
            pl.exclude(
                "replication", "cdr_rmsd", "cdr_rmsd_af2_full", "cdr_rmsd_af2_trim"
            )
        )
    )

    # # we must regenerate job name since this method changes the TCRs
    # val = serial_apply(
    #     val,
    #     shorten_both_tcrs,
    # )

    # val = generate_job_name(
    #     val,
    #     [
    #         "peptide",
    #         "mhc_1_seq",
    #         "mhc_2_seq",
    #         "tcr_1_seq",
    #         "tcr_2_seq",
    #     ],
    # )

    val.write_parquet(args.output_triad_path)

    pdb_antigen = val.select(FORMAT_ANTIGEN_COLS).unique()
    pdb_antigen = generate_job_name(
        pdb_antigen,
        ["peptide", "mhc_1_seq", "mhc_2_seq"],
    )

    pdb_antigen.select(["job_name"] + FORMAT_ANTIGEN_COLS).write_parquet(
        args.output_pmhc_path,
    )
