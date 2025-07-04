#!/usr/bin/env python3
"""
Script to calculate interface similarity between complexes and their MSA alignments
"""
import argparse
import polars as pl
from pathlib import Path
from tcrtrifold.interface_similarity import add_interface_similarity_features


def main():
    parser = argparse.ArgumentParser(
        description="Calculate interface similarity features"
    )
    parser.add_argument(
        "--data_dir", type=str, required=True, help="Base data directory"
    )
    parser.add_argument("--dset_name", type=str, required=True, help="Dataset name")
    parser.add_argument("--outdir", type=str, help="Output directory")

    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    dset_dir = data_dir / args.dset_name
    outdir = Path(args.outdir) if args.outdir else dset_dir / "interface_similarity"
    outdir.mkdir(parents=True, exist_ok=True)

    msa_dir = dset_dir / "msa"

    # Process p:MHC complexes
    pmhc_files = list((dset_dir / "pmhc" / "staged").glob("*.cleaned*.parquet"))
    for pmhc_file in pmhc_files:
        print(f"Processing p:MHC file: {pmhc_file}")
        df = pl.read_parquet(pmhc_file)

        df_with_features = add_interface_similarity_features(
            df, msa_dir, complex_type="pmhc"
        )

        output_file = outdir / f"{pmhc_file.stem}_interface_similarity.parquet"
        df_with_features.write_parquet(output_file)

        print(
            f"  - Processed {df.shape[0]} rows, added {len(df_with_features.columns) - len(df.columns)} interface similarity features"
        )
        print(f"  - Saved to: {output_file}")

    # Process triad complexes
    triad_files = list((dset_dir / "triad" / "staged").glob("*.cleaned*.parquet"))
    for triad_file in triad_files:
        print(f"Processing triad file: {triad_file}")
        df = pl.read_parquet(triad_file)

        # Process in smaller batches to avoid memory issues
        batch_size = 10
        all_results = []

        for i in range(0, len(df), batch_size):
            batch = df.slice(i, batch_size)
            batch_with_features = add_interface_similarity_features(
                batch, msa_dir, complex_type="triad"
            )
            all_results.append(batch_with_features)

        df_with_features = pl.concat(all_results)

        output_file = outdir / f"{triad_file.stem}_interface_similarity.parquet"
        df_with_features.write_parquet(output_file)

        print(
            f"  - Processed {df.shape[0]} rows, added {len(df_with_features.columns) - len(df.columns)} interface similarity features"
        )
        print(f"  - Saved to: {output_file}")


if __name__ == "__main__":
    main()
