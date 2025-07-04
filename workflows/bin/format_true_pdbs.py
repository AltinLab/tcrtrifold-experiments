import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pdb_parquet",
        type=str,
    )
    parser.add_argument(
        "--output_dir",
        type=str,
    )
    args = parser.parse_args()
