#!/usr/bin/env python

from tcrtrifold.utils import generate_job_name, FORMAT_ANTIGEN_COLS, TCRDIST_COLS
from tcrtrifold.tcr import shorten_tcr_to_vregion
import polars as pl
import argparse
import requests
from io import StringIO
import time
import re
from datetime import datetime, timezone
from Bio.Align import substitution_matrices
import Bio.Align

aligner = Bio.Align.PairwiseAligner(scoring="blastp")
BLOSUM = substitution_matrices.load("BLOSUM62")
AF3_CUTOFF = datetime(2023, 1, 12, tzinfo=timezone.utc)

PDB_QUERY = """
    query($id: String!) {
      entry(entry_id: $id) {
        polymer_entities {
          rcsb_polymer_entity_container_identifiers {
            entity_id
            asym_ids
            auth_asym_ids
          }
          entity_poly {
            pdbx_seq_one_letter_code_can
          }
        }
        rcsb_accession_info {
            initial_release_date
         }
      }
    }
    """


def blast_pident(seq_1, seq_2):
    aln = sorted(aligner.align(seq_1, seq_2))[0]
    aligned_q, aligned_s = aln.aligned
    a1, a2 = aln
    alnlen = sum(x != "-" for x in a1)
    ident = sum(x == y for x, y in zip(a1, a2))
    pident = 100.0 * ident / alnlen if alnlen else 0.0
    return pident, (a1, a2)


def get_chains_from_pdb(pdb_id):

    r = requests.post(
        "https://data.rcsb.org/graphql",
        json={"query": PDB_QUERY, "variables": {"id": pdb_id}},
        timeout=120,
    )
    r.raise_for_status()
    raw_pdb_dat = r.json()["data"]["entry"]

    pdb_dat = {
        "release_date": raw_pdb_dat["rcsb_accession_info"]["initial_release_date"],
        "chain_seqs": [
            polymer["entity_poly"]["pdbx_seq_one_letter_code_can"]
            for polymer in raw_pdb_dat["polymer_entities"]
        ],
    }

    return pdb_dat


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--query_triad_parquet",
        type=str,
    )
    parser.add_argument("--pmhc_parquet", type=str)
    parser.add_argument(
        "--blast_mhc_1_tsv",
        type=str,
    )
    parser.add_argument(
        "--blast_mhc_2_tsv",
        type=str,
    )
    parser.add_argument(
        "--blast_tcr_1_tsv",
        type=str,
    )
    parser.add_argument(
        "--blast_tcr_2_tsv",
        type=str,
    )
    parser.add_argument("--exclude_af3_training_only", action="store_true")
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

    mhc_1_ident = pl.read_csv(args.blast_mhc_1_tsv, separator="\t")

    mhc_2_ident = pl.read_csv(args.blast_mhc_2_tsv, separator="\t")

    tcr_1_ident = pl.read_csv(args.blast_tcr_1_tsv, separator="\t")

    tcr_2_ident = pl.read_csv(args.blast_tcr_2_tsv, separator="\t")

    strikeout_jn = []
    matched_triad = []
    pmhc_in_validation_jn = []
    matched_pmhc = []
    matched_peptide_aln = []

    pdb_cache = {}

    for row in query_triad.iter_rows(named=True):

        focal_mhc_1 = mhc_1_ident.filter(pl.col("qseqid") == row["mhc_1_hash"])
        focal_pdb_id = focal_mhc_1.select("pdb").unique()

        if focal_mhc_1.height == 0:
            continue

        if row["mhc_class"] == "II":
            # two conditions:
            # - found a match above thresh
            # - it's with a PDB ID we already matched
            focal_mhc_2 = mhc_2_ident.filter(
                pl.col("qseqid") == row["mhc_2_hash"]
            ).join(focal_pdb_id.select("pdb"), on="pdb")
            focal_pdb_id = focal_mhc_2.select("pdb").unique()

            if focal_mhc_2.height == 0:
                continue

        # now, parse the PDB IDs currently present
        pdb_matches = []
        peptide_matches = []
        for pdb_id in focal_pdb_id.select("pdb").to_series():
            if pdb_id not in pdb_cache:
                pdb_dat = get_chains_from_pdb(pdb_id)
                pdb_cache[pdb_id] = pdb_dat
            else:
                pdb_dat = pdb_cache[pdb_id]

            if args.exclude_af3_training_only:
                if (
                    datetime.strptime(
                        pdb_dat["release_date"], "%Y-%m-%dT%H:%M:%SZ"
                    ).replace(tzinfo=timezone.utc)
                    >= AF3_CUTOFF
                ):
                    continue

            seqs = pdb_dat["chain_seqs"]
            all_matches = []
            for seq in seqs:

                pident, aln_str_tuple = blast_pident(row["peptide"], seq)

                # 2 substitutions / 9 = 0.77
                if pident >= 77:
                    # heuristic to catch alignments to long chains that aren't occuring near
                    # the ends (bound by linker)
                    if len(seq) >= 27:
                        # capture first and last non-hyphen
                        first_non_hyph = (
                            re.compile(r"[^-]").search(aln_str_tuple[0])
                        ).regs[0][0]
                        last_non_hyph = (
                            re.compile(r"([^-])-*$").search(aln_str_tuple[0])
                        ).regs[0][0]

                        if not (first_non_hyph <= 10) and not (
                            last_non_hyph + len(row["peptide"])
                            >= len(aln_str_tuple[0]) - 10
                        ):
                            continue

                    all_matches.append((pident, aln_str_tuple))

            if len(all_matches) != 0:
                peptide_matches.append(all_matches[0][1][1])
                pdb_matches.append(pdb_id)

        if len(pdb_matches) != 0:
            pmhc_in_validation_jn.append(row["job_name"])
            matched_pmhc.append(pdb_matches)
            matched_peptide_aln.append(peptide_matches)

        focal_pdb_id = focal_pdb_id.filter(pl.col("pdb").is_in(pdb_matches))

        focal_tcr_1 = tcr_1_ident.filter(pl.col("qseqid") == row["tcr_1_hash"]).join(
            focal_pdb_id.select("pdb"), on="pdb"
        )

        if focal_tcr_1.height == 0:
            continue

        focal_pdb_id = focal_tcr_1.select("pdb").unique()

        focal_tcr_2 = tcr_2_ident.filter(pl.col("qseqid") == row["tcr_2_hash"]).join(
            focal_pdb_id.select("pdb"), on="pdb"
        )

        if focal_tcr_2.height == 0:
            continue

        strikeout_jn.append(row["job_name"])
        matched_triad.append(focal_tcr_2.select("pdb").to_series().to_list())

    pmhc_in_pdb = pl.DataFrame(
        {
            "job_name": pmhc_in_validation_jn,
            "pmhc_in_pdb": True,
            "pmhc_matches": matched_pmhc,
            "peptide_alignments": matched_peptide_aln,
        },
        schema={
            "job_name": pl.String,
            "pmhc_in_pdb": pl.Boolean,
            "pmhc_matches": pl.List(pl.String),
            "peptide_alignments": pl.List(pl.String),
        },
    ).unique()

    triad_in_pdb = pl.DataFrame(
        {
            "job_name": strikeout_jn,
            "triad_matches": matched_triad,
        },
        schema={"job_name": pl.String, "triad_matches": pl.List(pl.String)},
    ).unique()

    annot_triad = (
        query_triad.join(pmhc_in_pdb, on="job_name", how="left")
        .with_columns(
            pl.when(pl.col("pmhc_in_pdb").is_not_null())
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
            .alias("pmhc_in_pdb"),
            pl.when(pl.col("pmhc_in_pdb").is_not_null())
            .then(pl.col("pmhc_matches"))
            .otherwise(pl.lit(None))
            .alias("pmhc_matches"),
        )
        .select(
            pl.exclude(
                ["peptide_hash", "mhc_1_hash", "mhc_2_hash", "tcr_1_hash", "tcr_2_hash"]
            )
        )
    )

    exclude = query_triad.join(triad_in_pdb, on="job_name", how="inner").select(
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
