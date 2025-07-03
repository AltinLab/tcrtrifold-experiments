from .utils import (
    FORMAT_TCR_COLS,
    TCRDIST_COLS,
    FORMAT_ANTIGEN_COLS,
)
from tcrdist.repertoire import TCRrep
import polars as pl
import scipy


def pw_tcrdist(tcr_df, species=None, chains=None):

    if species is None:
        species = tcr_df[0].select("tcr_1_species").item()

    if chains is None:
        chain_1 = tcr_df[0].select("tcr_1_chain").item()
        chain_2 = tcr_df[0].select("tcr_2_chain").item()

        chains = [chain_1, chain_2]

    tcr_df_pd = tcr_df.with_columns(
        pl.col("tcr_1_cdr_1").alias("cdr1_a_aa"),
        pl.col("tcr_1_cdr_2").alias("cdr2_a_aa"),
        pl.col("tcr_1_cdr_2_5").alias("pmhc_a_aa"),
        pl.col("tcr_1_cdr_3").alias("cdr3_a_aa"),
        pl.col("tcr_2_cdr_1").alias("cdr1_b_aa"),
        pl.col("tcr_2_cdr_2").alias("cdr2_b_aa"),
        pl.col("tcr_2_cdr_2_5").alias("pmhc_b_aa"),
        pl.col("tcr_2_cdr_3").alias("cdr3_b_aa"),
        pl.lit(1).alias("count"),
    ).to_pandas()

    tr = TCRrep(
        cell_df=tcr_df_pd,
        organism=species,
        chains=chains,
        db_file="alphabeta_gammadelta_db.tsv",
        infer_cdrs=False,
    )

    # assert tr.pw_alpha.shape[0] == len(tcr_df_pd)

    dists = tr.pw_alpha + tr.pw_beta

    out_df = tr.clone_df.copy()

    out_df["index"] = out_df.index

    out_df.drop(
        labels=[
            "cdr3_a_aa",
            "cdr3_b_aa",
            "cdr1_a_aa",
            "cdr2_a_aa",
            "pmhc_a_aa",
            "cdr1_b_aa",
            "cdr2_b_aa",
            "pmhc_b_aa",
            "count",
            "clone_id",
        ],
        axis=1,
        inplace=True,
    )

    return pl.DataFrame(out_df), dists


def per_antigen_tcrdist_clust(df):

    df_by_antigen = df.partition_by(
        FORMAT_ANTIGEN_COLS,
    )

    out_dfs = []

    for antigen_df in df_by_antigen:

        tcr_df = antigen_df.select(FORMAT_TCR_COLS + TCRDIST_COLS)

        tcr_with_idx, pw_dist = pw_tcrdist(
            tcr_df,
            use_provided_cdr=True,
        )

        compressed = scipy.spatial.distance.squareform(pw_dist)
        Z = scipy.cluster.hierarchy.linkage(
            compressed,
            method="complete",
        )

        clusters = scipy.cluster.hierarchy.fcluster(
            Z,
            t=120,
            criterion="distance",
        )

        # add cluster labels to tcr_df
        tcr_with_idx = tcr_with_idx.with_columns(pl.Series("cluster", clusters))

        # add cluster labels to antigen_df
        antigen_df_clust = antigen_df.join(
            tcr_with_idx.select(FORMAT_TCR_COLS + TCRDIST_COLS + ["cluster"]),
            on=FORMAT_TCR_COLS + TCRDIST_COLS,
        )

        # now rank clusters by size, with largest cluster first
        cluster_sizes = (
            antigen_df_clust.group_by("cluster")
            .len()
            .sort("len", descending=False)
            .with_row_index(name="rank")
        )

        antigen_df_clust = antigen_df_clust.join(
            cluster_sizes,
            on="cluster",
        )

        out_dfs.append(antigen_df_clust)

    # combine all antigen dfs
    out_df = pl.concat(out_dfs)
    return out_df
