import polars as pl
from mdaf3.FeatureExtraction import serial_apply
import sqlite3
from tcrtrifold.utils import (
    generate_job_name,
    update_df_from_k_v,
    TCRDIST_COLS,
    FORMAT_COLS,
)
from tcrtrifold.cleaning import infer_hla_chain
from tcrtrifold.mhc import (
    B2M_HUMAN_SEQ,
    HLACodeWebConverter,
)
from tcrtrifold.tcr import extract_tcrdist_cols
from Bio import SeqIO
from io import StringIO
import requests

addtl_cols = TCRDIST_COLS + [
    "receptor_id",
    "references",
]

DB_PATH = "/tgen_labs/altin/alphafold3/IEDB/2025-04-15/iedb_public.db"

TRIAD_QUERY_STR = """
SELECT complex.*, 
    tcell.tcell_id, mhc_allele_name, mhc_allele_restriction.class as mhc_class,
    epitope.linear_peptide_seq
FROM complex
    JOIN reference on complex.reference_id = reference.reference_id
    JOIN tcell on tcell.complex_id = complex.complex_id
    JOIN mhc_allele_restriction on mhc_allele_restriction.mhc_allele_restriction_id = tcell.mhc_allele_restriction_id
    LEFT JOIN tcell_receptor ON tcell.tcell_id = tcell_receptor.tcell_id
    JOIN curated_epitope on tcell.curated_epitope_id = curated_epitope.curated_epitope_id
    JOIN epitope_object on curated_epitope.e_object_id = epitope_object.object_id
    JOIN epitope on epitope.epitope_id = epitope_object.epitope_id
WHERE tcell_receptor.tcell_id IS NULL;
"""

SEQ_STRUCT = pl.Struct(
    {
        "tcr_1_seq": pl.String,
        "tcr_2_seq": pl.String,
    }
)

pre_fasta_corrections = {
    "7Q9B": {
        "tcr_c1_pdb_chain": "D",
        "tcr_c2_pdb_chain": "E",
    },
    "6RPBA": {
        "pdb_id": "6RPB",
    },
}

# final_corrections = {
#     "6BJ3": {
#         "tcr_1_v_gene": "",
#         "tcr_2_j_gene": "",
#     },
#     "6BJ8": {
#         "tcr_1_v_gene": "",
#         "tcr_2_j_gene": "",
#     },
#     "8ES8": {
#         "tcr_1_v_gene": "",
#         "tcr_2_j_gene": "",
#     },
#     "8YIV": {
#         "tcr_1_v_gene": "",
#         "tcr_2_j_gene": "",
#     },
# }

invalid_pdb_ids = [
    # synthetic
    "9DL1",
    # in validation set
    "8ES9",
    "8I5C",
    "8EN8",
    # identified as delta alpha TCR
    # by anarci and igBlast
    "6BJ3",
    "6BJ8",
    "8ES8",
    "8YIV",
    "8YJ2",
]


def parse_chain(chain):
    if "[" in chain:

        return chain.split("[auth ")[1][0]
        # if can have multi-letter chains
        # return chain.split("[auth ")[1].split("]")[0]
    else:
        return chain.replace(" ", "")


def parse_fasta_description(description):
    chain_token = description.split("|")[1]

    if chain_token.startswith("Chain "):
        return list(parse_chain(chain_token.split("Chain ")[1]))
    else:
        chains = chain_token.split("Chains ")[1].split(",")
        chain_list = [parse_chain(chain) for chain in chains]

        return chain_list


def get_fasta_seq(
    pdb_id,
    tcr_1_chain_id,
    tcr_2_chain_id,
):
    r = requests.get("https://www.rcsb.org/fasta/entry/" + pdb_id)

    r.raise_for_status()

    fasta_sequences = SeqIO.parse(StringIO(r.text), "fasta")

    seq_dict = {}
    for fasta in fasta_sequences:
        chains = parse_fasta_description(fasta.description)
        for chain in chains:
            seq_dict[chain] = str(fasta.seq)

    try:

        return {
            "tcr_1_seq": seq_dict[tcr_1_chain_id] if tcr_1_chain_id is not None else "",
            "tcr_2_seq": seq_dict[tcr_2_chain_id] if tcr_2_chain_id is not None else "",
        }
    except KeyError:
        print("debug")


human_conv = HLACodeWebConverter()


conn = sqlite3.connect(DB_PATH)
iedb_metadat = (
    pl.read_database(
        connection=conn,
        query=TRIAD_QUERY_STR,
        infer_schema_length=1000000,
    )
    # .select("receptor_id", "assay_type")
    .unique()
)

iedb_metadat = iedb_metadat.filter(~pl.col("pdb_id").is_in(pl.Series(invalid_pdb_ids)))

iedb_I = (
    iedb_metadat.filter(pl.col("mhc_class") == "I")
    .rename({"mhc_allele_name": "mhc_1_name", "linear_peptide_seq": "peptide"})
    .with_columns(
        pl.col("mhc_1_name").str.split_exact("HLA-", 1).alias("split_parts"),
    )
    .select(pl.exclude("mhc_1_name"))
    .unnest("split_parts")
    .rename(
        {
            "field_0": "tmp",
            "field_1": "mhc_1_name",
        }
    )
    .select(pl.exclude("tmp", "mhc_struct", "split_parts", "mhc_name_inferred"))
    .with_columns(
        pl.lit("B2M").alias("mhc_2_name"),
        pl.lit(B2M_HUMAN_SEQ).alias("mhc_2_seq"),
        pl.lit("alpha").alias("mhc_1_chain"),
        pl.lit("beta").alias("mhc_2_chain"),
        pl.lit("human").alias("mhc_1_species"),
        pl.lit("human").alias("mhc_2_species"),
        pl.lit("alpha").alias("tcr_1_chain"),
        pl.lit("beta").alias("tcr_2_chain"),
        pl.lit("human").alias("tcr_1_species"),
        pl.lit("human").alias("tcr_2_species"),
    )
    .with_columns(
        pl.col("mhc_1_name")
        .map_elements(
            lambda x: human_conv.get_sequence(x, top_only=True),
            return_dtype=pl.String,
        )
        .alias("mhc_1_seq")
    )
)

for pdb_id, correction in pre_fasta_corrections.items():

    for k, v in correction.items():
        iedb_I = update_df_from_k_v(
            iedb_I,
            "pdb_id",
            pdb_id,
            k,
            v,
        )

iedb_I = iedb_I.with_columns(
    pl.struct(
        pl.col("pdb_id"),
        pl.col("tcr_c1_pdb_chain"),
        pl.col("tcr_c2_pdb_chain"),
    )
    .map_elements(
        lambda x: get_fasta_seq(
            x["pdb_id"],
            x["tcr_c1_pdb_chain"],
            x["tcr_c2_pdb_chain"],
        ),
        return_dtype=SEQ_STRUCT,
        skip_nulls=False,
    )
    .alias("chain_seqs"),
).unnest("chain_seqs")

iedb_II = (
    (
        iedb_metadat.filter(pl.col("mhc_class") == "II")
        .rename({"mhc_allele_name": "mhc_1_name", "linear_peptide_seq": "peptide"})
        .with_columns(pl.col("mhc_1_name").str.split("/").alias("split_parts"))
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
    )
    .with_columns(
        pl.lit("alpha").alias("mhc_1_chain"),
        pl.lit("beta").alias("mhc_2_chain"),
        pl.lit("alpha").alias("tcr_1_chain"),
        pl.lit("beta").alias("tcr_2_chain"),
        pl.lit("human").alias("mhc_1_species"),
        pl.lit("human").alias("mhc_2_species"),
        pl.lit("human").alias("tcr_1_species"),
        pl.lit("human").alias("tcr_2_species"),
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
    .with_columns(
        pl.struct(
            pl.col("pdb_id"),
            pl.col("tcr_c1_pdb_chain"),
            pl.col("tcr_c2_pdb_chain"),
        )
        .map_elements(
            lambda x: get_fasta_seq(
                x["pdb_id"],
                x["tcr_c1_pdb_chain"],
                x["tcr_c2_pdb_chain"],
            ),
            return_dtype=SEQ_STRUCT,
            skip_nulls=False,
        )
        .alias("chain_seqs"),
    )
    .unnest("chain_seqs")
    .select(pl.exclude("mhc_struct", "split_parts", "mhc_name_inferred"))
)

print(iedb_I.columns, iedb_II.columns)

iedb_metadat = pl.concat([iedb_I, iedb_II], how="align")

iedb_metadat = serial_apply(
    iedb_metadat,
    extract_tcrdist_cols,
)

iedb_metadat = generate_job_name(
    iedb_metadat,
    [
        "peptide",
        "mhc_1_seq",
        "mhc_2_seq",
        "tcr_1_seq",
        "tcr_2_seq",
    ],
)

iedb_metadat = iedb_metadat.with_columns(
    (
        pl.col("pdb_id").map_elements(
            lambda x: [[x]], return_dtype=pl.List(pl.List(pl.String))
        )
    ).alias("references"),
    pl.lit([["missing"]]).alias("receptor_id"),
    pl.lit(True).alias("cognate"),
)

iedb_metadat.select(
    FORMAT_COLS
    + addtl_cols
    + ["tcr_1_v_gene", "tcr_1_j_gene", "tcr_2_v_gene", "tcr_2_j_gene"]
).write_parquet("addtl_xray.parquet")
