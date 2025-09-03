import polars as pl
import hashlib
from Bio import SeqIO

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


def fasta_to_polars(fasta_path: str, desc_as_name: bool = False) -> pl.DataFrame:
    """
    Read a FASTA file and convert it into a Polars DataFrame
    with columns ["name", "sequence"].

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    Returns
    -------
    pl.DataFrame
        - "name": the sequence ID
        - "sequence": the full sequence string
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))

    if desc_as_name:
        names = [rec.description for rec in records]
    else:
        names = [rec.id for rec in records]
    seqs = [str(rec.seq) for rec in records]

    df = pl.DataFrame({"name": names, "seq": seqs})
    return df


def a3m_to_polars(a3m_handle, desc_as_name=False) -> pl.DataFrame:
    """
    Read an A3M (FASTA‐style) alignment from a file path or file‐like handle
    and convert it into a Polars DataFrame with columns ["name", "seq"].

    Parameters
    ----------
    a3m_handle : str or file‐like
        Path to the A3M file, or an open file handle / StringIO.
    desc_as_name : bool
        If True, use the full FASTA description line as the name;
        otherwise use only the record.id.

    Returns
    -------
    pl.DataFrame
        - "name": the sequence identifier (or full description)
        - "seq":  the raw sequence string (including lowercase insertions)
    """
    # If given a file path, open it; otherwise assume handle semantics
    if isinstance(a3m_handle, str):
        handle = open(a3m_handle, "r")
        close_when_done = True
    else:
        handle = a3m_handle
        close_when_done = False

    # Parse all records
    records = list(SeqIO.parse(handle, "fasta"))

    if close_when_done:
        handle.close()

    # Extract names and sequences
    if desc_as_name:
        names = [rec.description for rec in records]
    else:
        names = [rec.id for rec in records]
    seqs = [str(rec.seq) for rec in records]

    # Build and return DataFrame
    return pl.DataFrame({"name": names, "seq": seqs})


def filter_to_cog_thresh(triad_df, thresh, lt=False):
    if lt:
        antigen_count = (
            triad_df.filter(pl.col("cognate"))
            .group_by(FORMAT_ANTIGEN_COLS)
            .len(name="n_tcr")
            .filter(pl.col("n_tcr") <= thresh)
            .select(FORMAT_ANTIGEN_COLS + ["n_tcr"])
        )
    else:
        antigen_count = (
            triad_df.filter(pl.col("cognate"))
            .group_by(FORMAT_ANTIGEN_COLS)
            .len(name="n_tcr")
            .filter(pl.col("n_tcr") >= thresh)
            .select(FORMAT_ANTIGEN_COLS + ["n_tcr"])
        )
    triad_df = triad_df.join(antigen_count, on=FORMAT_ANTIGEN_COLS, how="inner")
    return triad_df
