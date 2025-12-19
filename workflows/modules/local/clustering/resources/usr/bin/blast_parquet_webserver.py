#!/usr/bin/env python

from tcrtrifold.utils import generate_job_name, FORMAT_ANTIGEN_COLS, TCRDIST_COLS
from tcrtrifold.tcr import shorten_tcr_to_vregion
import polars as pl
import argparse
import requests
from io import StringIO
import time
import re

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

RCSB_API = "https://search.rcsb.org/rcsbsearch/v2/query"


def blast_rcsb(seq_id, seq):
    json = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "sequence_type": "protein",
                "value": seq,
                "identity_cutoff": 0.95,
            },
        },
        "request_options": {"paginate": {"start": 0, "rows": 10000}},
        "return_type": "entry",
    }

    r = requests.post(RCSB_API, json=json, timeout=10000)
    r.raise_for_status()

    if r.status_code == 204:
        return pl.DataFrame(
            schema={"pdb": pl.String, "score": pl.Float64, "qseqid": pl.String}
        )

    data_dict = r.json()["result_set"]

    return (
        pl.DataFrame(data_dict)
        .with_columns(pl.lit(seq_id).alias("qseqid"))
        .rename({"identifier": "pdb"})
    )


# deprecated: NCBI blast against PDB misses obvious hits
def blastp_pdb_df(seq, evalue=1e-3, max_hits=1000, poll_s=3):
    # Submit job
    r = requests.post(
        BLAST_URL,
        data={
            "CMD": "Put",
            "PROGRAM": "blastp",
            "DATABASE": "pdb",
            "QUERY": seq,
            "EXPECT": evalue,
            "HITLIST_SIZE": max_hits,
        },
    )
    r.raise_for_status()

    rid = re.compile(r"RID = (.+)").search(r.text).group(1)

    # Poll until ready
    while True:
        r = requests.get(
            BLAST_URL, params={"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"}
        )
        r.raise_for_status()

        status = re.compile(r"Status=(.+)").search(r.text).group(1)
        if status == "READY":
            break
        if "FAILED" == status or "UNKNOWN" == status:
            raise RuntimeError(f"BLAST job {rid} failed or expired")
        time.sleep(poll_s)

    # Get tabular results (outfmt 6)
    res = requests.get(
        BLAST_URL,
        params={
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": "CSV",
            "ALIGNMENT_VIEW": "Tabular",
            "DESCRIPTIONS": max_hits,
            "ALIGNMENTS": max_hits,
        },
    ).text

    # Load into Polars
    # https://rnnh.github.io/bioinfo-notebook/docs/blast.html outfmt
    schema = {
        "qseqid": pl.String,
        "sseqid": pl.String,
        "pident": pl.Float64,
        "length": pl.Int32,
        "mismatch": pl.Int32,
        "gapopen": pl.Int32,
        "qstart": pl.Int32,
        "qend": pl.Int32,
        "sstart": pl.Int32,
        "send": pl.Int32,
        "evalue": pl.Float64,
        "bitscore": pl.Float64,
        # "percent positive"
        "ppos": pl.Float64,
    }
    col_names = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "ppos",
    ]
    df = pl.read_csv(
        StringIO(res), has_header=False, new_columns=col_names, schema=schema
    ).filter(pl.col("qseqid").is_not_null())
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--seq_parquet",
        type=str,
    )
    parser.add_argument("--output_blast_tsv")
    args = parser.parse_args()

    # # read fasta into string
    # with open(args.fasta, "r") as f:
    #     fasta_content = f.read()

    # # blast_df = blastp_pdb_df(
    # #     fasta_content,
    # # )

    seq_parquet = pl.read_parquet(args.seq_parquet)

    seq_dfs = []
    for row in seq_parquet.iter_rows(named=True):
        seq_dfs.append(blast_rcsb(row["name"], row["seq"]))

    pl.concat(seq_dfs).write_csv(args.output_blast_tsv, separator="\t")
