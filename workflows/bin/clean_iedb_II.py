#!/usr/bin/env python
from tcrtrifold.utils import generate_job_name
from tcrtrifold.cleaning import infer_hla_chain
from tcrtrifold.mhc import (
    B2M_HUMAN_SEQ,
    HLACodeWebConverter,
)
from tcrtrifold.utils import (
    generate_job_name,
    FORMAT_COLS,
    FORMAT_TCR_COLS,
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

    human_conv = HLACodeWebConverter()

    addtl_cols = TCRDIST_COLS + [
        "receptor_id",
        "references",
    ]

    schema_overrides = {
        "references": pl.String,
        "receptor_id": pl.String,
    }

    iedb_human_II = pl.read_csv(args.raw_csv_path, schema_overrides=schema_overrides)

    iedb_human_II = (
        iedb_human_II.rename(
            {
                "Peptide": "peptide",
                "TCRb": "tcr_2_seq",
                "TCRa": "tcr_1_seq",
                "HLA": "mhc_1_name",
            }
        )
        .with_columns(
            pl.lit("alpha").alias("mhc_1_chain"),
            pl.lit("beta").alias("mhc_2_chain"),
            pl.lit("alpha").alias("tcr_1_chain"),
            pl.lit("beta").alias("tcr_2_chain"),
            pl.lit("human").alias("tcr_1_species"),
            pl.lit("human").alias("tcr_2_species"),
            pl.lit("human").alias("mhc_1_species"),
            pl.lit("human").alias("mhc_2_species"),
            pl.lit("II").alias("mhc_class"),
            pl.lit(True).alias("cognate"),
            pl.col("receptor_id").str.split(",").alias("receptor_id"),
            pl.col("references").str.split(",").alias("references"),
        )
        .filter(
            pl.col("tcr_1_seq").is_not_null(),
            pl.col("tcr_2_seq").is_not_null(),
        )
        .with_columns(
            pl.when(pl.col("references").is_not_null())
            .then(
                pl.col("references").list.eval(pl.element().str.split("/").list.get(-1))
            )
            .otherwise(None)
            .alias("references")
        )
    )

    iedb_human_II = (
        iedb_human_II.with_columns(
            pl.col("mhc_1_name").str.split("/").alias("split_parts")
        )
        .with_columns(
            pl.when(pl.col("split_parts").list.len() == 2)
            .then(
                pl.struct(
                    pl.col("split_parts")
                    .list.get(0, null_on_oob=True)
                    .str.slice(4)
                    .alias("mhc_1_name"),
                    pl.col("split_parts")
                    .list.get(1, null_on_oob=True)
                    .alias("mhc_2_name"),
                )
            )
            .otherwise(
                pl.struct(
                    pl.lit(None).alias("mhc_1_name"),
                    pl.col("split_parts").list.get(0).str.slice(4).alias("mhc_2_name"),
                )
            )
            .alias("mhc_struct")
        )
        .select(pl.exclude("mhc_1_name"))
        .with_columns(
            pl.col("mhc_struct")
            .map_elements(
                lambda x: infer_hla_chain(x["mhc_1_name"], x["mhc_2_name"]),
                return_dtype=pl.Struct(
                    {
                        "mhc_1_name": pl.String,
                        "mhc_2_name": pl.String,
                        "mhc_name_inferred": pl.String,
                    }
                ),
            )
            .alias("chains")
        )
        .unnest("chains")
        .filter(
            (pl.col("mhc_1_name").is_not_null()) & (pl.col("mhc_2_name").is_not_null())
        )
        .with_columns(
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
    )

    iedb_human_II = generate_job_name(
        iedb_human_II,
        [
            "peptide",
            "mhc_1_seq",
            "mhc_2_seq",
            "tcr_1_seq",
            "tcr_2_seq",
        ],
    )

    iedb_human_II = serial_apply(
        iedb_human_II,
        extract_tcrdist_cols,
    )

    iedb_human_II = iedb_human_II.group_by(FORMAT_COLS + TCRDIST_COLS).agg(
        [pl.col("references").flatten(), pl.col("receptor_id").flatten()]
    )

    iedb_human_II.select(FORMAT_COLS + addtl_cols).write_parquet(
        args.output_triad_path,
    )

    iedb_human_II_antigen = iedb_human_II.select(FORMAT_ANTIGEN_COLS).unique()
    iedb_human_II_antigen = generate_job_name(
        iedb_human_II_antigen,
        ["peptide", "mhc_1_seq", "mhc_2_seq"],
    )

    iedb_human_II_antigen.select(["job_name"] + FORMAT_ANTIGEN_COLS).write_parquet(
        args.output_pmhc_path,
    )
