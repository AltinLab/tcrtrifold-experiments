#!/usr/bin/env python
from tcrtrifold.mhc import (
    HLASequenceDBConverter,
    H2SequenceDictConverter,
)
from tcrtrifold.utils import (
    generate_job_name,
    update_df_from_k_v,
    FORMAT_COLS,
    FORMAT_TCR_COLS,
    FORMAT_ANTIGEN_COLS,
    TCRDIST_COLS,
)
from tcrtrifold.tcr import (
    extract_tcrdist_cols,
)
from mdaf3.FeatureExtraction import serial_apply, split_apply_combine
import polars as pl
import argparse
from datetime import datetime, timezone
import requests
import warnings
from Bio import SeqIO
from io import StringIO

SEQ_STRUCT = pl.Struct(
    {
        "peptide": pl.String,
        "mhc_1_seq": pl.String,
        "mhc_2_seq": pl.String,
        "tcr_1_seq": pl.String,
        "tcr_2_seq": pl.String,
    }
)


def download_pdb(row, path):

    r = requests.get(f"https://files.rcsb.org/download/{row["pdb"]}.pdb")
    suffix = ".pdb"
    try:
        r.raise_for_status()
    except Exception as e:
        r = requests.get(f"https://files.rcsb.org/download/{row["pdb"]}.cif")
        suffix = ".cif"
        r.raise_for_status()
    with open(path / (row["pdb"] + suffix), "wb") as f:
        f.write(r.content)
    return row


def extract_pdb_date(row):
    r = requests.get("https://data.rcsb.org/rest/v1/core/entry/" + row["pdb"])
    r.raise_for_status()
    new_row = row.copy()
    new_row["pdb_date"] = r.json()["rcsb_accession_info"]["initial_release_date"]

    return new_row


def get_pdb_date(df):
    df = split_apply_combine(
        df,
        extract_pdb_date,
        chunksize=50,
    ).with_columns(pl.col("pdb_date").str.to_datetime().alias("pdb_date"))
    return df


# these serve as a log of changes we make to the PDB data in order to get it into our format
pre_fasta_corrections = {
    "3tf7": {
        "mhc_1_segid": "E",
        "tcr_1_segid": "G",
        "tcr_2_segid": "G",
        "peptide_segid": "F",
    },
    # antigen chain located on MHC 2 / TCR 2 chain
    # when we later look for the peptide we can find it by aligning to this chain
    "6bga": {"peptide_segid": "B"},
    "3pl6": {"peptide_segid": "D"},
    "3o6f": {"peptide_segid": "B"},
    "6dfw": {"peptide_segid": "D"},
    "3c5z": {"peptide_segid": "H"},
    "3rdt": {"peptide_segid": "D"},
    "6dfx": {"peptide_segid": "B"},
    "3c60": {"peptide_segid": "D"},
    "4grl": {"peptide_segid": "D"},
    "6mnn": {"peptide_segid": "D"},
    "6dfs": {"peptide_segid": "D"},
    "4p4k": {"peptide_segid": "B"},
    "4may": {"peptide_segid": "D"},
    "8vd0": {"peptide_segid": "C"},
    "4qrp": {
        "mhc_1_segid": "A",
        "mhc_2_segid": "B",
        "peptide_segid": "C",
        "tcr_1_segid": "D",
        "tcr_2_segid": "E",
    },
}

post_fasta_corrections = {
    # 3tf7 has 1 tcrab pair bound together with linker
    "3tf7": {
        "tcr_1_seq": "MGAQSVTQPDARVTVSEGASLQLRCKYSYSATPYLFWYVQYPRQGPQMLLKYYSGDPVVQGVNGFEAEFSKSDSSFHLRKASVHRSDSAVYFCAVSAKGTGSKLSFGKGAKLTVSP",
        "tcr_2_seq": "SEAAVTQSPRNKVTVTGENVTLSCRQTNSHNYMYWYRQDTGHELRLIYYSYGAGNLQIGDVPDGYKATRTTQEDFFLTLESASPSQTSLYFCASSDAPGQLYFGEGSKLTVLELEHHHHHH",
    },
    "8vd0": {"peptide": "GQVELGGGNAVEVCKG"},
    "7q9b": {
        "mhc_1_segid": "FFF",
        "mhc_2_segid": "GGG",
        "tcr_1_segid": "III",
        "tcr_2_segid": "JJJ",
        "peptide_segid": "HHH",
        "mhc_1_species": "human",
    },
    "6bga": {"peptide": "ADSLSFFSSSIKR"},
    "3pl6": {"peptide": "MKENPVVHFFKNIVTPR"},
    "3o6f": {"peptide": "FSWGAEGQRPGFGS"},
    "6dfw": {"peptide": "HLVERLYLVCGGEGAGGG"},
    "3c5z": {"peptide": "FEAQKAKANKAVD"},
    "3rdt": {"peptide": "FEAQKAKANKAVD"},
    "6dfx": {"peptide": "VEELYLVAGEEGCG"},
    "3c60": {"peptide": "FEAQKAKANKAVD"},
    "4grl": {"peptide": "MKDRLLMLFAKDVVSRN"},
    "6mnn": {"peptide": "RVSYYGPKTSPVQ"},
    "6dfs": {"peptide": "HLVERLYLVCGEEGA"},
    "4p4k": {"peptide": "QAFWIDLFETIGGGSLVPRGS"},
    "4may": {"peptide": "MKRQLVHFVRDFAQLGGS"},
}

# this is a record of PDB IDs post AF3 cutoff that we find acceptably formatted for our use
post_af3_cutoff_valid_pdbs = [
    "8gom",
    "8vd0",
    "8trr",
    "8wte",
    "8vcy",
    "8es9",
    "8gon",
    "8i5d",
    "8i5c",
    "8vcx",
    "7q99",
    "8eo8",
    "8dnt",
    "8enh",
    "8ye4",
    "8wul",
    "8f5a",
    "8vd2",
    "7q9b",
    "8en8",
    "8pjg",
    "7q9a",
]

# this PDB ID was not found by STCR, but was included in Phil Bradley's 2023 paper
outlier_pdb = pl.DataFrame(
    [
        {
            "pdb": "6l9l",
            "mhc_1_segid": "A",
            "mhc_2_segid": None,
            "peptide_segid": "B",
            "tcr_1_segid": "C",
            "tcr_2_segid": "D",
            "mhc_1_species": "mouse",
            "mhc_2_species": None,
            "tcr_1_species": "mouse",
            "tcr_2_species": "mouse",
            "mhc_class": "I",
            "mhc_1_chain": "heavy",
            "mhc_2_chain": None,
            "cognate": True,
            "tcr_1_chain": "alpha",
            "tcr_2_chain": "beta",
            "pdb_date": extract_pdb_date({"pdb": "6l9l"})["pdb_date"],
        }
    ]
).with_columns(pl.col("pdb_date").str.to_datetime().alias("pdb_date"))


def format_stcr_df(df):
    """
    STCR dat dataframes contain duplicate PDB rows
    Aggregate them and fail out on unexpected cases
    """

    # probably some mislabeling, but "GA"/"GB" and "MHC2" are used interchangeably in "mhc_type"
    # column
    df = df.with_columns(
        pl.when((pl.col("mhc_type") == "GA") | (pl.col("mhc_type") == "GB"))
        .then(pl.lit("MH2"))
        .otherwise(pl.col("mhc_type"))
        .alias("mhc_type"),
    )

    # sometimes, antigen_type is null. In all cases we have seen, this is a peptide antigen
    df = df.with_columns(
        pl.when(pl.col("antigen_type").is_null())
        .then(pl.lit("peptide"))
        .otherwise(pl.col("antigen_type"))
        .alias("antigen_type"),
    )

    # sometimes stcr finds multiple antigens per triad, some non-peptide
    # find the segid of the peptide antigen in the "|" separated list of peptide
    # segids (orig colname antigen_chain)
    df = (
        df.filter(pl.col("antigen_type").str.contains("peptide"))
        .with_columns(
            pl.col("antigen_type")
            .str.split("|")
            .list.eval(pl.element().str.strip_chars())
            .list.eval(pl.element().index_of("peptide"))
            .list.first()
            .alias("peptide_segid_index")
        )
        .with_columns(
            pl.col("antigen_chain")
            .str.split("|")
            .list.eval(pl.element().str.strip_chars())
            .list.get(pl.col("peptide_segid_index"))
            .alias("antigen_chain")
        )
    )

    # df = df.with_columns(
    #     pl.concat_list(
    #         [
    #             "mhc_chain1",
    #             "mhc_chain2",
    #             "antigen_chain",
    #             "Achain",
    #             "Bchain",
    #         ]
    #     ).alias("segid_list")
    # )

    df = df.group_by(pl.col("pdb"), maintain_order=True).agg(
        pl.col("mhc_type").drop_nulls(),
        pl.col("mhc_chain1").first().alias("mhc_1_segid"),
        pl.col("mhc_chain2").first().alias("mhc_2_segid"),
        pl.col("antigen_chain").first().alias("peptide_segid"),
        pl.col("Achain").first().alias("tcr_1_segid"),
        pl.col("Bchain").first().alias("tcr_2_segid"),
        # pl.col("segid_list"),
        pl.col("mhc_chain1_organism").drop_nulls().alias("mhc_1_species"),
        pl.col("mhc_chain2_organism").drop_nulls().alias("mhc_2_species"),
        pl.col("alpha_organism").drop_nulls().alias("tcr_1_species"),
        pl.col("beta_organism").drop_nulls().alias("tcr_2_species"),
    )

    if df.filter(pl.col("mhc_type").list.n_unique() > 1).height > 0:
        raise ValueError(
            "Unexpected MHC type ambiguity in STCR data. "
            "Please check the input data for inconsistencies."
        )

    # reformat mhc_type -> mhc_class
    df = df.with_columns(
        pl.col("mhc_type").list.first().alias("mhc_type"),
    )

    df = df.with_columns(
        pl.when(pl.col("mhc_type") == "MH1")
        .then(pl.lit("I"))
        .when(pl.col("mhc_type") == "MH2")
        .then(pl.lit("II"))
        .otherwise(None)
        .alias("mhc_class"),
    ).drop("mhc_type")

    # organism error handling

    org_cols = [
        "mhc_1_species",
        "mhc_2_species",
        "tcr_1_species",
        "tcr_2_species",
    ]

    for org_col in org_cols:
        if df.filter(pl.col(org_col).list.n_unique() > 1).height > 0:
            raise ValueError(
                f"Unexpected species ambiguity in {org_col}. "
                "Please check the input data for inconsistencies"
            )

        df = df.with_columns(pl.col(org_col).list.first().alias(org_col))

    # df = df.filter(
    #     (pl.col("mhc_chain1").is_not_null())
    #     & (pl.col("mhc_chain2").is_not_null())
    # )

    # df = df.group_by("pdb").agg(
    #     pl.col("Bchain").drop_nulls().first(),
    #     pl.col("Achain").drop_nulls().first(),
    #     pl.col("mhc_chain1").drop_nulls().first(),
    #     pl.col("mhc_chain2").drop_nulls().first(),
    #     pl.col("antigen_chain").drop_nulls().first(),
    #     pl.col("mhc_class").drop_nulls().first(),
    #     pl.col("mhc_chain1_organism").drop_nulls().first().alias("mhc_1_species"),
    #     pl.col("mhc_chain2_organism").drop_nulls().first().alias("mhc_2_species"),
    #     pl.col("alpha_organism").drop_nulls().first().alias("tcr_1_species"),
    #     pl.col("beta_organism").drop_nulls().first().alias("tcr_2_species"),
    # )

    df = df.with_columns(
        pl.when(pl.col("mhc_1_species") == "homo sapiens")
        .then(pl.lit("human"))
        .when(pl.col("mhc_1_species") == "mus musculus")
        .then(pl.lit("mouse"))
        .otherwise(None)
        .alias("mhc_1_species"),
        pl.when(pl.col("mhc_2_species") == "homo sapiens")
        .then(pl.lit("human"))
        .when(pl.col("mhc_2_species") == "mus musculus")
        .then(pl.lit("mouse"))
        .otherwise(None)
        .alias("mhc_2_species"),
        pl.when(pl.col("tcr_1_species") == "homo sapiens")
        .then(pl.lit("human"))
        .when(pl.col("tcr_1_species") == "mus musculus")
        .then(pl.lit("mouse"))
        .otherwise(None)
        .alias("tcr_1_species"),
        pl.when(pl.col("tcr_2_species") == "homo sapiens")
        .then(pl.lit("human"))
        .when(pl.col("tcr_2_species") == "mus musculus")
        .then(pl.lit("mouse"))
        .otherwise(None)
        .alias("tcr_2_species"),
    )

    df = df.with_columns(
        pl.when(pl.col("mhc_class") == "II")
        .then(pl.lit("alpha"))
        .otherwise(pl.lit("heavy"))
        .alias("mhc_1_chain"),
        pl.when(pl.col("mhc_class") == "II")
        .then(pl.lit("beta"))
        .otherwise(pl.lit("light"))
        .alias("mhc_2_chain"),
        pl.lit(True).alias("cognate"),
        pl.lit("alpha").alias("tcr_1_chain"),
        pl.lit("beta").alias("tcr_2_chain"),
    )

    return df


def refmt_rep_dset(rep):
    rep = rep.rename({"pdbid": "pdb"}).with_columns(
        pl.when(pl.col("mhc_class") == 1)
        .then(pl.lit("I"))
        .otherwise(pl.lit("II"))
        .alias("mhc_class"),
    )

    return rep


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
    antigen_chain_id,
    mhc_chain1_id,
    mhc_chain2_id,
    Achain_id,
    Bchain_id,
):
    r = requests.get("https://www.rcsb.org/fasta/entry/" + pdb_id)

    r.raise_for_status()

    fasta_sequences = SeqIO.parse(StringIO(r.text), "fasta")

    seq_dict = {}
    for fasta in fasta_sequences:
        chains = parse_fasta_description(fasta.description)
        for chain in chains:
            seq_dict[chain] = str(fasta.seq)

    return {
        "peptide": (seq_dict[antigen_chain_id] if antigen_chain_id is not None else ""),
        "mhc_1_seq": (seq_dict[mhc_chain1_id] if mhc_chain1_id is not None else ""),
        "mhc_2_seq": (seq_dict[mhc_chain2_id] if mhc_chain2_id is not None else ""),
        "tcr_1_seq": seq_dict[Achain_id] if Achain_id is not None else "",
        "tcr_2_seq": seq_dict[Bchain_id] if Bchain_id is not None else "",
    }


def format_seqs(df, skip_peptide=False):
    df = (
        df.with_columns(
            pl.struct(
                pl.col("pdb"),
                pl.col("mhc_1_segid"),
                pl.col("mhc_2_segid"),
                pl.col("peptide_segid"),
                pl.col("tcr_1_segid"),
                pl.col("tcr_2_segid"),
            )
            .map_elements(
                lambda x: get_fasta_seq(
                    x["pdb"],
                    x["peptide_segid"],
                    x["mhc_1_segid"],
                    x["mhc_2_segid"],
                    x["tcr_1_segid"],
                    x["tcr_2_segid"],
                ),
                return_dtype=SEQ_STRUCT,
                skip_nulls=False,
            )
            .alias("chain_seqs"),
        )
        .unnest("chain_seqs")
        .with_columns(
            pl.when(pl.col("peptide") == "")
            .then(pl.lit(None))
            .otherwise(pl.col("peptide"))
            .alias("peptide"),
            pl.when(pl.col("mhc_1_seq") == "")
            .then(pl.lit(None))
            .otherwise(pl.col("mhc_1_seq"))
            .alias("mhc_1_seq"),
            pl.when(pl.col("mhc_2_seq") == "")
            .then(pl.lit(None))
            .otherwise(pl.col("mhc_2_seq"))
            .alias("mhc_2_seq"),
            pl.when(pl.col("tcr_1_seq") == "")
            .then(pl.lit(None))
            .otherwise(pl.col("tcr_1_seq"))
            .alias("tcr_1_seq"),
            pl.when(pl.col("tcr_2_seq") == "")
            .then(pl.lit(None))
            .otherwise(pl.col("tcr_2_seq"))
            .alias("tcr_2_seq"),
        )
    )

    return df


def remove_peptide_from_chains(row):
    new_row = row.copy()

    if row["mhc_1_seq"] is not None and row["peptide"] in row["mhc_1_seq"]:
        warnings.warn(
            f"Peptide found in MHC 1 sequence for PDB {row['pdb']} at position {row['mhc_1_seq'].index(row['peptide'])}"
        )
        index_of_peptide = row["mhc_1_seq"].index(row["peptide"])
        new_row["mhc_1_seq"] = new_row["mhc_1_seq"][
            index_of_peptide + len(row["peptide"]) :
        ]
    if row["mhc_2_seq"] is not None and row["peptide"] in row["mhc_2_seq"]:
        warnings.warn(
            f"Peptide found in MHC 2 sequence for PDB {row['pdb']} at position {row['mhc_2_seq'].index(row['peptide'])}"
        )
        index_of_peptide = row["mhc_2_seq"].index(row["peptide"])
        new_row["mhc_2_seq"] = new_row["mhc_2_seq"][
            index_of_peptide + len(row["peptide"]) :
        ]
    if row["tcr_1_seq"] is not None and row["peptide"] in row["tcr_1_seq"]:
        warnings.warn(
            f"Peptide found in TCR 1 sequence for PDB {row['pdb']} at position {row['tcr_1_seq'].index(row['peptide'])}"
        )
        index_of_peptide = row["tcr_1_seq"].index(row["peptide"])
        new_row["tcr_1_seq"] = new_row["tcr_1_seq"][
            index_of_peptide + len(row["peptide"]) :
        ]
    if row["tcr_2_seq"] is not None and row["peptide"] in row["tcr_2_seq"]:
        warnings.warn(
            f"Peptide found in TCR 2 sequence for PDB {row['pdb']} at position {row['tcr_2_seq'].index(row['peptide'])}"
        )
        index_of_peptide = row["tcr_2_seq"].index(row["peptide"])
        new_row["tcr_2_seq"] = new_row["tcr_2_seq"][
            index_of_peptide + len(row["peptide"]) :
        ]
    return new_row


def infer_correct_mhc(row, human_conv, mouse_conv):
    mhc1 = row["mhc_1_seq"]
    mhc2 = row["mhc_2_seq"]

    if row["mhc_1_species"] == "human":
        mhc_1_inf = human_conv.get_mhc_allele(
            mhc1, chain=row["mhc_1_chain"], top_only=False
        )
    else:
        mhc_1_inf = mouse_conv.get_mhc_allele(
            mhc1, chain=row["mhc_1_chain"], top_only=False
        )

    if row["mhc_1_species"] == "human":
        if row["mhc_class"] == "I" and row["mhc_2_seq"] is None:
            mhc_2_inf = {
                "seq": None,
                "name": None,
                "match_size": None,
                "max_resolution_name": None,
                "sequence_status": None,
            }
        else:
            mhc_2_inf = human_conv.get_mhc_allele(
                mhc2, chain=row["mhc_2_chain"], top_only=False
            )
    else:
        if row["mhc_class"] == "I" and row["mhc_2_seq"] is None:
            mhc_2_inf = {
                "seq": None,
                "name": None,
                "match_size": None,
                "max_resolution_name": None,
                "sequence_status": None,
            }
        else:
            mhc_2_inf = mouse_conv.get_mhc_allele(
                mhc2, chain=row["mhc_2_chain"], top_only=False
            )

    new_row = row.copy()

    new_row["mhc_1_match_seq"] = mhc_1_inf["seq"]
    new_row["mhc_1_name"] = mhc_1_inf["name"]
    new_row["mhc_1_match_size"] = mhc_1_inf["match_size"]
    new_row["mhc_1_match_proportion"] = (
        (mhc_1_inf["match_size"] / len(mhc1))
        if mhc_1_inf["match_size"] is not None
        else None
    )
    new_row["mhc_1_status"] = mhc_1_inf["sequence_status"]
    new_row["mhc_1_name_maxres"] = mhc_1_inf["max_resolution_name"]

    new_row["mhc_2_match_seq"] = mhc_2_inf["seq"]
    new_row["mhc_2_name"] = mhc_2_inf["name"]
    new_row["mhc_2_match_size"] = mhc_2_inf["match_size"]
    new_row["mhc_2_match_proportion"] = (
        (mhc_2_inf["match_size"] / len(mhc2))
        if mhc_2_inf["match_size"] is not None
        else None
    )
    new_row["mhc_2_status"] = mhc_2_inf["sequence_status"]
    new_row["mhc_2_maxres"] = mhc_2_inf["max_resolution_name"]
    return new_row


# def refmt_stcr_dest(stcr):
#     stcr = format_pdb_df(stcr)
#     stcr = get_pdb_date(stcr)

#     return stcr


def refmt_mhc_name(mhc_name_df):
    mhc_name_df_II = (
        mhc_name_df.filter(pl.col("mhc_class") == "II")
        .with_columns(pl.col("mhc").str.split(",").alias("split_parts"))
        .with_columns(
            pl.when(pl.col("split_parts").list.len() == 2)
            .then(
                pl.struct(
                    pl.col("split_parts")
                    .list.get(0, null_on_oob=True)
                    .alias("mhc_1_name"),
                    pl.col("split_parts")
                    .list.get(1, null_on_oob=True)
                    .alias("mhc_2_name"),
                )
            )
            .otherwise(
                pl.struct(
                    pl.lit(None).alias("mhc_1_name"),
                    pl.col("split_parts").list.get(0).alias("mhc_2_name"),
                )
            )
            .alias("mhc_struct")
        )
        .unnest("mhc_struct")
    )

    mhc_name_df_I = (
        mhc_name_df.filter(pl.col("mhc_class") == "I")
        .with_columns(pl.col("mhc").str.split(",").alias("split_parts"))
        .with_columns(
            pl.when(pl.col("split_parts").list.len() == 2)
            .then(
                pl.struct(
                    pl.col("split_parts")
                    .list.get(0, null_on_oob=True)
                    .alias("mhc_1_name"),
                    pl.col("split_parts")
                    .list.get(1, null_on_oob=True)
                    .alias("mhc_2_name"),
                )
            )
            .otherwise(
                pl.struct(
                    pl.lit("B2M").alias("mhc_2_name"),
                    pl.col("split_parts").list.get(0).alias("mhc_1_name"),
                )
            )
            .alias("mhc_struct")
        )
        .unnest("mhc_struct")
    )

    return pl.concat(
        [
            mhc_name_df_II,
            mhc_name_df_I,
        ],
        how="vertical",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--raw_csv_path",
        type=str,
    )
    parser.add_argument(
        "-d",
        "--raw_stcr_path",
        type=str,
    )
    parser.add_argument(
        "-i",
        "--imgt_hla_path",
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

    # https://www.nature.com/articles/s41586-024-07487-w#data-availability
    cutoff = pl.lit(datetime(2023, 1, 12, tzinfo=timezone.utc))

    schema_overrides = {
        "Gchain": pl.String,
        "Dchain": pl.String,
    }
    null_values = ["NA", "unknown", "NOT"]

    stcr_raw = pl.read_csv(
        args.raw_stcr_path,
        schema_overrides=schema_overrides,
        null_values=null_values,
        separator="\t",
    )

    stcr = format_stcr_df(stcr_raw)
    stcr = get_pdb_date(stcr)

    rep = pl.read_csv(
        args.raw_csv_path,
    )

    rep = refmt_rep_dset(rep)

    # we keep a PDB ID for two reasons:
    # 1. it is in Phil Bradley's 2023 paper
    # 2. it is in PDB post-AF3 cutoff and is useable by our pipeline

    rep_pdb_ls = rep.select("pdb").to_series().to_list()

    keep_pdb = rep_pdb_ls + post_af3_cutoff_valid_pdbs

    keep_df = stcr.filter(pl.col("pdb").is_in(keep_pdb))

    keep_df = pl.concat([keep_df, outlier_pdb])

    keep_df = keep_df.with_columns(
        pl.when(pl.col("pdb").is_in(post_af3_cutoff_valid_pdbs))
        .then(pl.lit(False))
        .otherwise(pl.lit(True))
        .alias("replication")
    )

    # keep_df = keep_df.join(
    #     rep.select(
    #         [
    #             "pdb",
    #             "organism",
    #         ]
    #     ),
    #     on="pdb",
    #     how="left",
    # )

    # keep_df = keep_df.with_columns(
    #     pl.when(pl.col("organism").is_null())
    #     .then(pl.lit("human"))
    #     .otherwise(pl.col("organism"))
    #     .alias("tcr_1_species"),
    #     pl.when(pl.col("organism").is_null())
    #     .then(pl.lit("human"))
    #     .otherwise(pl.col("organism"))
    #     .alias("tcr_2_species"),
    # )

    # use the organism information from Phil Bradley's paper
    # the new PDBs happen to all be human

    # keep_df = keep_df.with_columns(
    #     pl.when(pl.col("organism").is_null())
    #     .then(pl.lit("human"))
    #     .otherwise(pl.col("organism"))
    #     .alias("organism"),
    # )

    for pdb_id, correction in pre_fasta_corrections.items():

        for k, v in correction.items():
            keep_df = update_df_from_k_v(
                keep_df,
                "pdb",
                pdb_id,
                k,
                v,
            )

    # now extract seqs from FASTA
    keep_df = format_seqs(keep_df)

    # in some cases, peptides were manually extracted in Phil Bradley's paper
    # we keep these for consistency with the paper
    # keep_df = (
    #     keep_df.join(
    #         rep.select(
    #             "pdb",
    #             "peptide",
    #         ).rename({"peptide": "peptide_manually_extracted"}),
    #         how="left",
    #         on="pdb",
    #     )
    #     .with_columns(
    #         pl.when(pl.col("peptide_manually_extracted").is_not_null())
    #         .then(pl.col("peptide_manually_extracted"))
    #         .otherwise(pl.col("peptide"))
    #         .alias("peptide")
    #     )
    #     .drop(
    #         "peptide_manually_extracted",
    #     )
    # )

    for pdb_id, correction in post_fasta_corrections.items():
        for k, v in correction.items():
            keep_df = update_df_from_k_v(
                keep_df,
                "pdb",
                pdb_id,
                k,
                v,
            )

    keep_df = keep_df.with_columns(
        pl.when(pl.col("tcr_1_species").is_null())
        .then(pl.col("mhc_1_species"))
        .otherwise(pl.col("tcr_1_species"))
        .alias("tcr_1_species"),
        pl.when(pl.col("tcr_2_species").is_null())
        .then(pl.col("mhc_1_species"))
        .otherwise(pl.col("tcr_2_species"))
        .alias("tcr_2_species"),
    )

    # rep_mhc_names = refmt_mhc_name(rep.select("pdb", "mhc_class", "mhc"))

    # new_mhc_names = serial_apply(
    #     keep_df.filter(~pl.col("replication")),
    #     infer_correct_mhc,
    #     HLASequenceDBConverter(args.imgt_hla_path),
    #     None,
    # ).rename(
    #     {
    #         "mhc_1_name": "mhc_1_name_inferred",
    #         "mhc_2_name": "mhc_2_name_inferred",
    #     }
    # )

    keep_df = serial_apply(
        keep_df,
        infer_correct_mhc,
        HLASequenceDBConverter(args.imgt_hla_path),
        H2SequenceDictConverter(),
    )

    # keep_df = keep_df.join(
    #     new_mhc_names.select(
    #         [
    #             "pdb",
    #             "mhc_1_name_inferred",
    #             "mhc_2_name_inferred",
    #         ]
    #     ),
    #     on="pdb",
    #     how="left",
    # )

    # keep_df = keep_df.with_columns(
    #     pl.when(pl.col("mhc_1_name").is_null())
    #     .then(pl.col("mhc_1_name_inferred"))
    #     .otherwise(pl.col("mhc_1_name"))
    #     .alias("mhc_1_name"),
    #     pl.when(pl.col("mhc_2_name").is_null())
    #     .then(pl.col("mhc_2_name_inferred"))
    #     .otherwise(pl.col("mhc_2_name"))
    #     .alias("mhc_2_name"),
    # )

    keep_df = keep_df.join(
        rep.select(
            [
                "pdb",
                "cdr_rmsd",
                "cdr_rmsd_af2_full",
                "cdr_rmsd_af2_trim",
            ]
        ),
        on="pdb",
        how="left",
    )

    keep_df = generate_job_name(
        keep_df,
        [
            "peptide",
            "mhc_1_seq",
            "mhc_2_seq",
            "tcr_1_seq",
            "tcr_2_seq",
        ],
    )

    keep_df = serial_apply(
        keep_df,
        extract_tcrdist_cols,
    )

    keep_df = serial_apply(keep_df, remove_peptide_from_chains)

    keep_df.select(
        FORMAT_COLS
        + TCRDIST_COLS
        + [
            "pdb",
            "pdb_date",
            "replication",
            "peptide_segid",
            "mhc_1_segid",
            "mhc_2_segid",
            "tcr_1_segid",
            "tcr_2_segid",
            "cdr_rmsd",
            "cdr_rmsd_af2_full",
            "cdr_rmsd_af2_trim",
        ]
    ).write_parquet(
        args.output_triad_path,
    )

    pdb_antigen = keep_df.select(FORMAT_ANTIGEN_COLS).unique()
    pdb_antigen = generate_job_name(
        pdb_antigen,
        ["peptide", "mhc_1_seq", "mhc_2_seq"],
    )

    pdb_antigen.select(["job_name"] + FORMAT_ANTIGEN_COLS).write_parquet(
        args.output_pmhc_path,
    )
