import polars as pl
import hashlib

FORMAT_COLS = [
    "job_name",
    "cognate",
    "peptide",
    "mhc_class",
    "mhc_1_chain",
    "mhc_1_species",
    "mhc_1_name",
    "mhc_1_seq",
    "mhc_2_chain",
    "mhc_2_species",
    "mhc_2_name",
    "mhc_2_seq",
    "tcr_1_chain",
    "tcr_1_species",
    "tcr_1_seq",
    "tcr_2_chain",
    "tcr_2_species",
    "tcr_2_seq",
]

FORMAT_TCR_COLS = [
    "tcr_1_chain",
    "tcr_1_species",
    "tcr_1_seq",
    "tcr_2_chain",
    "tcr_2_species",
    "tcr_2_seq",
]

FORMAT_ANTIGEN_COLS = [
    "peptide",
    "mhc_class",
    "mhc_1_chain",
    "mhc_1_species",
    "mhc_1_name",
    "mhc_1_seq",
    "mhc_2_chain",
    "mhc_2_species",
    "mhc_2_name",
    "mhc_2_seq",
]

FORMAT_MHC_COLS = [
    "mhc_class",
    "mhc_1_chain",
    "mhc_1_species",
    "mhc_1_name",
    "mhc_1_seq",
    "mhc_2_chain",
    "mhc_2_species",
    "mhc_2_name",
    "mhc_2_seq",
]

TCRDIST_COLS = [
    "tcr_1_cdr_1",
    "tcr_1_cdr_2",
    "tcr_1_cdr_2_5",
    "tcr_1_cdr_3",
    "tcr_2_cdr_1",
    "tcr_2_cdr_2",
    "tcr_2_cdr_2_5",
    "tcr_2_cdr_3",
]


SOURCE_ANTIGEN_COLS = ["source_" + colname for colname in FORMAT_ANTIGEN_COLS]
SOURCE_RENAME_DICT = {k: v for k, v in zip(FORMAT_ANTIGEN_COLS, SOURCE_ANTIGEN_COLS)}
SOURCE_REV_RENAME_DICT = {
    v: k for k, v in zip(FORMAT_ANTIGEN_COLS, SOURCE_ANTIGEN_COLS)
}


def hash_sequence(seq: str, hash_type: str = "md5") -> str:
    """
    Hash a TCR sequence using the specified hash function.

    Args:
        tcr_seq (str): The TCR sequence string.
        hash_type (str): The hash function to use ('md5', 'sha1', 'sha256', etc.)

    Returns:
        str: The hexadecimal digest of the hashed sequence.
    """
    # Select the hash function
    if hash_type.lower() == "md5":
        h = hashlib.md5()
    elif hash_type.lower() == "sha1":
        h = hashlib.sha1()
    elif hash_type.lower() == "sha256":
        h = hashlib.sha256()
    else:
        raise ValueError("Unsupported hash type")

    # Encode the sequence and compute the hash
    h.update(seq.encode("utf-8"))
    return h.hexdigest()


def generate_job_name(df, cols, name="job_name"):
    df = df.with_columns(
        pl.concat_str(
            pl.concat_str(
                [
                    *[pl.col(colname) for colname in cols],
                ],
                ignore_nulls=True,
            )
            .map_elements(lambda x: hash_sequence(x, "md5"), return_dtype=pl.String)
            .alias(name),
        )
    )
    return df


def update_df_from_k_v(
    df,
    primary_key_colname,
    primary_key,
    k,
    v,
):
    df = pl.concat(
        [
            df.filter(pl.col(primary_key_colname) == primary_key).with_columns(
                pl.lit(v).alias(k)
            ),
            df.filter(pl.col(primary_key_colname) != primary_key),
        ],
        how="vertical_relaxed",
    )
    return df
