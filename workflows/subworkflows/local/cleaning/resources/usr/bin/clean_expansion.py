#!/usr/bin/env python
"""
stitchrdl -s human


"""
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
    FORMAT_MHC_COLS,
)
from tcrtrifold.mhc import DQA_FOR, DPA_FOR, DRA_NAME, DRA_SEQ, HLACodeWebConverter
from tcrtrifold.tcr import (
    shorten_tcr_to_vregion,
    extract_tcrdist_cols,
)
import tidytcells as tt

from mdaf3.FeatureExtraction import serial_apply
import polars as pl
import argparse
import warnings

# from Stitchr import stitchrfunctions as fxn
# from Stitchr import stitchr as st

# # Stitchr warnings
# warnings.filterwarnings(
#     "error",
#     message="Note: while a C-terminal CDR3:J germline match has been found, it was only the string*",
# )


# # stitchr initialization
# # https://jamieheather.github.io/stitchr/importing.html
# HUMAN_ALPHA_TCRDAT, HUMAN_ALPHA_FUNC, HUMAN_ALPHA_PARTIAL = fxn.get_ref_data(
#     "TRA", st.gene_types, "HUMAN"
# )
# HUMAN_BETA_TCRDAT, HUMAN_BETA_FUNC, HUMAN_BETA_PARTIAL = fxn.get_ref_data(
#     "TRB", st.gene_types, "HUMAN"
# )
# HUMAN_CODONS = fxn.get_optimal_codons("", "HUMAN")
# HUMAN_J_RES, LOW_CONF_JS = fxn.get_j_motifs("HUMAN")
# HUMAN_C_RES = fxn.get_c_motifs("HUMAN")


# def stitch_tcr(chain, v_gene, j_gene, cdr_3):
#     tcr_bits = {
#         "v": v_gene,
#         "j": j_gene,
#         "cdr3": cdr_3,
#         "l": v_gene,
#         # infer
#         "c": "",
#         "species": "HUMAN",
#         "seamless": False,
#         "5_prime_seq": "",
#         "3_prime_seq": "",
#         "name": "TCR",
#         "no_leader": True,
#         "skip_n_checks": False,
#         "skip_c_checks": False,
#         "mode": "",
#     }
#     tcr_bits, _ = fxn.sort_input(tcr_bits)
#     try:
#         if chain == "alpha":
#             out = st.stitch(
#                 tcr_bits,
#                 HUMAN_ALPHA_TCRDAT,
#                 HUMAN_ALPHA_FUNC,
#                 HUMAN_ALPHA_PARTIAL,
#                 HUMAN_CODONS,
#                 3,
#                 "",
#                 HUMAN_C_RES,
#                 HUMAN_J_RES,
#                 LOW_CONF_JS,
#             )
#         else:
#             out = st.stitch(
#                 tcr_bits,
#                 HUMAN_BETA_TCRDAT,
#                 HUMAN_BETA_FUNC,
#                 HUMAN_BETA_PARTIAL,
#                 HUMAN_CODONS,
#                 3,
#                 "",
#                 HUMAN_C_RES,
#                 HUMAN_J_RES,
#                 LOW_CONF_JS,
#             )
#     except Exception as e:
#         print(f"Stitchr failed: {e}")
#         return None

#     return fxn.translate_nt(out["stitched_nt"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--clonotype_shortlist",
        type=str,
    )
    parser.add_argument("--tcr_info", type=str)
    parser.add_argument(
        "--hla_genotyping",
        type=str,
    )

    parser.add_argument(
        "--peptide_pool",
        type=str,
    )

    parser.add_argument(
        "--output_triad_path",
        type=str,
    )
    args = parser.parse_args()

    conv = HLACodeWebConverter()

    ctypes = pl.read_csv(args.clonotype_shortlist)
    tcr_info = pl.read_csv(args.tcr_info)

    tcr_info = (
        tcr_info.with_columns(
            pl.concat_str(
                ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
            ).alias("tcr_seq")
        )
        .select("tcr_seq", "cdr3_nt", "v_gene", "j_gene", "chain")
        .unique(maintain_order=True)
        .group_by("cdr3_nt", "v_gene", "j_gene", "chain", maintain_order=True)
        .agg(pl.col("tcr_seq").first())
    )

    hla_types = pl.read_csv(args.hla_genotyping, null_values=["X"])
    pep_pool = pl.read_csv(args.peptide_pool, null_values=["#N/A"])

    dpa_cols = [
        "HLA-DPA1.1",
        "HLA-DPA1.2",
        "HLA-DPA2.1",
        "HLA-DPA2.2",
    ]

    dpb_cols = [
        "HLA-DPB1.1",
        "HLA-DPB1.2",
        "HLA-DPB2.1",
        "HLA-DPB2.2",
    ]

    dqa_cols = [
        "HLA-DQA1.1",
        "HLA-DQA1.2",
    ]

    dqb_cols = [
        "HLA-DQB1.1",
        "HLA-DQB1.2",
    ]

    dra_cols = [
        "HLA-DRA.1",
        "HLA-DRA.2",
    ]

    drb_cols = [
        "HLA-DRB1.1",
        "HLA-DRB1.2",
        "HLA-DRB2.1",
        "HLA-DRB2.2",
        "HLA-DRB3.1",
        "HLA-DRB3.2",
        "HLA-DRB4.1",
        "HLA-DRB4.2",
        "HLA-DRB5.1",
        "HLA-DRB5.2",
        "HLA-DRB6.1",
        "HLA-DRB6.2",
        "HLA-DRB7.1",
        "HLA-DRB7.2",
        "HLA-DRB8.1",
        "HLA-DRB8.2",
        "HLA-DRB9.1",
        "HLA-DRB9.2",
    ]

    b_chain_cols = dpb_cols + dqb_cols + drb_cols

    focal_hla_type = hla_types.filter(
        pl.col("SampleID (provided on PBMC tube)") == "75064_E"
    )

    ctypes = ctypes.filter(
        pl.col("TRA").is_not_null()
        & pl.col("TRB").is_not_null()
        & ~pl.col("TRA").str.contains(";")
        & ~pl.col("TRB").str.contains(";")
    )

    ctypes = (
        ctypes.with_columns(
            pl.col("TRA").str.split(by=":"), pl.col("TRB").str.split(by=":")
        )
        .with_columns(
            pl.col("TRA").list.to_struct(
                fields=["tcr_1_v_gene", "tcr_1_cdr_3", "tcr_1_j_gene", "tcr_1_cdr3_nt"]
            ),
            pl.col("TRB").list.to_struct(
                fields=["tcr_2_v_gene", "tcr_2_cdr_3", "tcr_2_j_gene", "tcr_2_cdr3_nt"]
            ),
        )
        .unnest("TRA", "TRB")
    )

    ctypes = (
        ctypes.join(
            tcr_info.filter(pl.col("chain") == "TRA"),
            left_on=["tcr_1_cdr3_nt", "tcr_1_v_gene", "tcr_1_j_gene"],
            right_on=["cdr3_nt", "v_gene", "j_gene"],
        )
        .rename({"tcr_seq": "tcr_1_seq"})
        .join(
            tcr_info.filter(pl.col("chain") == "TRB"),
            left_on=["tcr_2_cdr3_nt", "tcr_2_v_gene", "tcr_2_j_gene"],
            right_on=["cdr3_nt", "v_gene", "j_gene"],
        )
        .rename({"tcr_seq": "tcr_2_seq"})
    )

    ctypes = (
        ctypes.with_columns(
            pl.lit("human").alias("tcr_1_species"),
            pl.lit("alpha").alias("tcr_1_chain"),
            pl.lit("human").alias("tcr_2_species"),
            pl.lit("beta").alias("tcr_2_chain"),
        )
        .select(FORMAT_TCR_COLS)
        .unique()
    )

    ctypes = serial_apply(
        ctypes,
        extract_tcrdist_cols,
    )

    a_chains = []
    a_chain_seq = []
    b_chains = []
    b_chain_seq = []

    for b_chain_col in b_chain_cols:
        b_chain_type = b_chain_col[4:8] + "*"

        b_chain_allele = focal_hla_type.select(b_chain_col).item()

        if b_chain_allele is None:
            continue

        b_chain_allele = b_chain_type + b_chain_allele[:5]

        if b_chain_type.startswith("DQB"):
            a_chain_name = DQA_FOR.get(b_chain_allele)

        elif b_chain_type.startswith("DPB"):
            a_chain_name = DPA_FOR.get(b_chain_allele)

        else:
            a_chain_name = DRA_NAME

        if a_chain_name is None:
            print(f"No known alpha chain for {b_chain_allele}")
        else:
            a_chains.append(a_chain_name)
            a_chain_seq.append(conv.get_sequence(a_chain_name, top_only=True))
            b_chains.append(b_chain_allele)
            b_chain_seq.append(conv.get_sequence(b_chain_allele, top_only=True))

    mhc_df = pl.DataFrame(
        {
            "mhc_class": ["II"] * len(a_chains),
            "mhc_1_chain": ["alpha"] * len(a_chains),
            "mhc_1_name": a_chains,
            "mhc_1_seq": a_chain_seq,
            "mhc_1_species": ["human"] * len(a_chains),
            "mhc_2_name": b_chains,
            "mhc_2_seq": b_chain_seq,
            "mhc_2_chain": ["beta"] * len(b_chains),
            "mhc_2_species": ["human"] * len(b_chains),
        }
    )

    triad_df = ctypes.join(
        mhc_df,
        how="cross",
    ).join(pep_pool.select("peptide").unique(), how="cross")

    triad_df = generate_job_name(
        triad_df,
        [
            "peptide",
            "mhc_1_seq",
            "mhc_2_seq",
            "tcr_1_seq",
            "tcr_2_seq",
        ],
    )

    triad_df = triad_df.with_columns(pl.lit(None).alias("cognate"))

    triad_df = triad_df.select(FORMAT_COLS + TCRDIST_COLS)

    triad_df.write_parquet(args.output_triad_path)

    triad_df.select()

    # for row in ctypes.iter_rows(named=True):

    # tcr_1_v_gene = tt.tr.standardize(row["tcr_1_v_gene"])
    # tcr_1_j_gene = tt.tr.standardize(row["tcr_1_j_gene"])
    # tcr_2_v_gene = tt.tr.standardize(row["tcr_2_v_gene"])
    # tcr_2_j_gene = tt.tr.standardize(row["tcr_2_j_gene"])

    # tcr_1_seq = stitch_tcr("alpha", tcr_1_v_gene, tcr_1_j_gene, row["tcr_1_cdr_3"])
    # tcr_2_seq = stitch_tcr("beta", tcr_2_v_gene, tcr_2_j_gene, row["tcr_2_cdr_3"])

    # if tcr_1_seq is None or tcr_2_seq is None:
    #     print("Failed to stitch")
    #     continue

    # tcr_1_seq = shorten_tcr_to_vregion(tcr_1_seq, "alpha", "human")
    # tcr_2_seq = shorten_tcr_to_vregion(tcr_2_seq, "beta", "human")
    # print(tcr_1_seq)
    # print(tcr_2_seq)
