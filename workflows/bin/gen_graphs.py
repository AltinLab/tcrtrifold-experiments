#!/usr/bin/env python
from tcrtrifold.graph_utils import mda_triad_to_graph
from mdaf3.AF3OutputParser import AF3Output
from mdaf3.BoltzOutputParser import BoltzOutput
from torch_geometric.data import OnDiskDataset
from pathlib import Path
import argparse
import polars as pl


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_parquet",
        type=str,
    )
    parser.add_argument("--inference_type")
    parser.add_argument("--inference_dir")
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
    )

    args = parser.parse_args()

    if args.inference_type == "af3":
        inf_cls = AF3Output
    elif args.inference_type == "boltz":
        inf_cls = BoltzOutput
    else:
        raise ValueError

    parent_inf = Path(args.inference_dir)

    triad_df = pl.read_parquet(args.input_parquet)

    graphs = []

    for row in triad_df.iter_rows(named=True):

        inf_handle = inf_cls(parent_inf / row["job_name"])

        graph = mda_triad_to_graph(
            inf_handle.get_mda_universe(),
            inf_handle.get_contact_prob_ndarr(),
            row["mhc_class"],
            row["cognate"],
            row["job_name"],
        )

        graphs.append(graph)

    out_dset = OnDiskDataset(args.output_path)
    out_dset.extend(graphs)
    out_dset.close()
