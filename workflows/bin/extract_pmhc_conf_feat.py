#!/usr/bin/env python
from tcrtrifold.feat_extract import (
    extract_mean_tcr_pmhc_pae,
    extract_num_contacts,
    extract_mean_tcr_pmhc_pae_class_II,
    extract_summary_metrics,
    extract_peptide_pLDDT,
    extract_peptide_pLDDT_class_II,
    extract_cdr_pLDDT,
    extract_mhc_helix_pLDDT,
    extract_mean_peptide_mhc_pae,
)
from mdaf3.FeatureExtraction import split_apply_combine
import polars as pl
from pathlib import Path
import argparse
import warnings

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_parquet",
        type=str,
    )
    parser.add_argument("--inference_type")
    parser.add_argument("--inference_dir")
    parser.add_argument("--summary_only", action="store_true", default=False)
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
    )

    args = parser.parse_args()

    inference_type = args.inference_type
    inf_dir = Path(args.inference_dir)

    df = pl.read_parquet(args.input_parquet)

    df = split_apply_combine(
        df,
        extract_summary_metrics,
        inf_dir,
        inference_type,
        chunksize=15,
    )

    if not args.summary_only:

        df = split_apply_combine(
            df,
            extract_peptide_pLDDT,
            inf_dir,
            inference_type,
            chunksize=15,
        )

        df = split_apply_combine(
            df,
            extract_mhc_helix_pLDDT,
            inf_dir,
            inference_type,
            chunksize=15,
        )

        df = split_apply_combine(
            df,
            extract_mean_peptide_mhc_pae,
            inf_dir,
            inference_type,
            chunksize=15,
        )
        # class-II specific features
        if df.select("mhc_class")[0].item() == "II":

            df = split_apply_combine(
                df,
                extract_peptide_pLDDT_class_II,
                inf_dir,
                inference_type,
                chunksize=15,
            )

    df.write_parquet(args.output_path)
