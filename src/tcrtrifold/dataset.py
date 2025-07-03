from .utils import FORMAT_ANTIGEN_COLS
import polars as pl


class TriadDataset:

    def __init__(self, triad_df, antigen_df):

        if isinstance(triad_df, str):
            triad_df = pl.read_parquet(triad_df)

        if isinstance(antigen_df, str):
            antigen_df = pl.read_parquet(antigen_df)

        self.triad_df = triad_df
        self.antigen_df = antigen_df


def filt_to_cognate_thresh(triad_dataset, thresh):
    triad = triad_dataset.triad_df
    antigen = triad_dataset.antigen_df

    valid_antigen = (
        (
            triad.filter(pl.col("cognate"))
            .group_by(FORMAT_ANTIGEN_COLS)
            .len()
            .filter(pl.col("len") >= thresh)
        )
        .select(FORMAT_ANTIGEN_COLS)
        .unique()
    )

    triad = triad.join(
        valid_antigen,
        on=FORMAT_ANTIGEN_COLS,
        how="inner",
    )

    antigen = antigen.join(
        valid_antigen,
        on=FORMAT_ANTIGEN_COLS,
        how="inner",
    )

    return TriadDataset(triad, antigen)
