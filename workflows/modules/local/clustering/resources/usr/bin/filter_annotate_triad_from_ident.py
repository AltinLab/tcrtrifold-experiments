#!/usr/bin/env python

from tcrtrifold.utils import generate_job_name, FORMAT_ANTIGEN_COLS, TCRDIST_COLS
import polars as pl
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--query_triad_parquet",
        type=str,
    )
    parser.add_argument("--pmhc_parquet", type=str)
    parser.add_argument(
        "--validation_triad_parquet",
        type=str,
    )
    parser.add_argument(
        "--peptide_ident_tsv",
    )
    parser.add_argument(
        "--mhc_1_ident_tsv",
    )
    parser.add_argument(
        "--mhc_2_ident_tsv",
    )
    parser.add_argument(
        "--tcr_1_ident_tsv",
    )
    parser.add_argument(
        "--tcr_2_ident_tsv",
    )
    parser.add_argument(
        "--output_annotated_triad_parquet",
        type=str,
    )
    parser.add_argument(
        "--output_excluded_triad_parquet",
        type=str,
    )
    parser.add_argument(
        "--output_remaining_pmhc_parquet",
        type=str,
    )
    args = parser.parse_args()

    query_triad = pl.read_parquet(args.query_triad_parquet)
    curr_pmhc = pl.read_parquet(args.pmhc_parquet)

    validation_triad = pl.read_parquet(args.validation_triad_parquet)

    ident_cols = "query,target,pident,alnlen,qcov,tcov,evalue,bits,qstart,qend,tstart,tend".split(
        ","
    )

    peptide_ident = pl.read_csv(
        args.peptide_ident_tsv, separator="\t", has_header=False, new_columns=ident_cols
    )
    mhc_1_ident = pl.read_csv(
        args.mhc_1_ident_tsv, separator="\t", has_header=False, new_columns=ident_cols
    )
    mhc_2_ident = pl.read_csv(
        args.mhc_2_ident_tsv, separator="\t", has_header=False, new_columns=ident_cols
    )
    tcr_1_ident = pl.read_csv(
        args.tcr_1_ident_tsv, separator="\t", has_header=False, new_columns=ident_cols
    )
    tcr_2_ident = pl.read_csv(
        args.tcr_2_ident_tsv, separator="\t", has_header=False, new_columns=ident_cols
    )

    strikeout_jn = []
    pmhc_in_validation_jn = []

    for row in query_triad.iter_rows(named=True):

        v_subset = validation_triad

        focal_pep = peptide_ident.filter(pl.col("query") == row["peptide_hash"])

        if focal_pep.height == 0:
            continue

        v_subset = v_subset.join(
            focal_pep.select("target", "pident").rename(
                {"target": "peptide_hash", "pident": "peptide_ident"}
            ),
            on="peptide_hash",
        )

        focal_mhc_1 = mhc_1_ident.filter(pl.col("query") == row["mhc_1_hash"])

        if focal_mhc_1.height == 0:
            continue

        v_subset = v_subset.join(
            focal_mhc_1.select("target", "pident").rename(
                {"target": "mhc_1_hash", "pident": "mhc_1_ident"}
            ),
            on="mhc_1_hash",
        )

        if row["mhc_class"] == "II":
            focal_mhc_2 = mhc_2_ident.filter(pl.col("query") == row["mhc_2_hash"])

            if focal_mhc_2.height == 0:
                continue

            v_subset = v_subset.join(
                focal_mhc_2.select("target", "pident").rename(
                    {"target": "mhc_2_hash", "pident": "mhc_2_ident"}
                ),
                on="mhc_2_hash",
            )

            if (
                v_subset.filter(
                    pl.col("peptide_ident") >= 77,
                    pl.col("mhc_1_ident") >= 95,
                    pl.col("mhc_2_ident") >= 95,
                ).height
                > 0
            ):
                pmhc_in_validation_jn.append(row["job_name"])

        else:

            if (
                v_subset.filter(
                    pl.col("peptide_ident") >= 77, pl.col("mhc_1_ident") >= 95
                ).height
                > 0
            ):
                pmhc_in_validation_jn.append(row["job_name"])

        focal_tcr_1 = tcr_1_ident.filter(pl.col("query") == row["tcr_1_hash"])

        if focal_tcr_1.height == 0:
            continue

        v_subset = v_subset.join(
            focal_tcr_1.select("target", "pident").rename(
                {"target": "tcr_1_hash", "pident": "tcr_1_ident"}
            ),
            on="tcr_1_hash",
        )

        focal_tcr_2 = tcr_2_ident.filter(pl.col("query") == row["tcr_2_hash"])

        if focal_tcr_2.height == 0:
            continue

        v_subset = v_subset.join(
            focal_tcr_2.select("target", "pident").rename(
                {"target": "tcr_2_hash", "pident": "tcr_2_ident"}
            ),
            on="tcr_2_hash",
        )

        if v_subset.height != 0:

            if row["mhc_class"] == "II":
                if (
                    v_subset.filter(
                        pl.col("peptide_ident") >= 77,
                        pl.col("mhc_1_ident") >= 95,
                        pl.col("mhc_2_ident") >= 95,
                        pl.col("tcr_1_ident") >= 95,
                        pl.col("tcr_2_ident") >= 95,
                    ).height
                    > 0
                ):
                    strikeout_jn.append(row["job_name"])
            else:

                if (
                    v_subset.filter(
                        pl.col("peptide_ident") >= 80,
                        pl.col("mhc_1_ident") >= 95,
                        pl.col("tcr_1_ident") >= 95,
                        pl.col("tcr_2_ident") >= 95,
                    ).height
                    > 0
                ):
                    strikeout_jn.append(row["job_name"])

    pmhc_in_v = pl.DataFrame(
        {
            "job_name": pmhc_in_validation_jn,
            "pmhc_in_validation": True,
        }
    ).unique()

    triad_in_v = pl.DataFrame(
        {
            "job_name": strikeout_jn,
        },
        schema={"job_name": pl.String},
    ).unique()

    annot_triad = (
        query_triad.join(pmhc_in_v, on="job_name", how="left")
        .with_columns(
            pl.when(pl.col("pmhc_in_validation").is_not_null())
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias("pmhc_in_validation"),
        )
        .select(
            pl.exclude(
                ["peptide_hash", "mhc_1_hash", "mhc_2_hash", "tcr_1_hash", "tcr_2_hash"]
            )
        )
    )

    exclude = query_triad.join(triad_in_v, on="job_name", how="inner").select(
        pl.exclude(
            [
                "peptide_hash",
                "mhc_1_hash",
                "mhc_2_hash",
                "tcr_1_hash",
                "tcr_2_hash",
            ]
        )
    )

    annot_triad = annot_triad.join(
        exclude.select("job_name").unique(), on="job_name", how="anti"
    )

    remaining_antigen = curr_pmhc.join(
        generate_job_name(
            annot_triad.select(FORMAT_ANTIGEN_COLS).unique(),
            ["peptide", "mhc_1_seq", "mhc_2_seq"],
            name="job_name",
        )
        .select("job_name")
        .unique(),
        on="job_name",
    )

    annot_triad.write_parquet(
        args.output_annotated_triad_parquet,
    )

    exclude.write_parquet(
        args.output_excluded_triad_parquet,
    )
    remaining_antigen.write_parquet(
        args.output_remaining_pmhc_parquet,
    )
