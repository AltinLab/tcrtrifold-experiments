from .utils import (
    FORMAT_COLS,
    FORMAT_TCR_COLS,
    FORMAT_ANTIGEN_COLS,
    FORMAT_MHC_COLS,
    FORMAT_MHC_COLS,
    TCRDIST_COLS,
    SOURCE_ANTIGEN_COLS,
    SOURCE_RENAME_DICT,
    SOURCE_REV_RENAME_DICT,
    generate_job_name,
)
from .tcrdist_utils import pw_tcrdist
from tqdm import tqdm
import editdistance
import numpy as np
import polars as pl
import random


def generate_all_possible_negs(
    pos_df,
    negs_from_df,
    tcrdist_thresh=120,
    pmhc_edit_thresh=3,
):
    """
    Generate all possible negatives from `negs_from_df` for the antigens
    in `pos_df`.

    Both DFs must have cols: FORMAT_COLS + TCRDIST_COLS

    A sampled negative TCR for an antigen in `pos_df` is only valid if:
    - it comes from a negative antigen (peptide + mhc) >= `pmhc_edit_thresh` edit distance
        from the `pos_df` antigen
    - it, and all the TCRs that come from the negative antigen, are >= `tcrdist_thresh`
        from all the TCRs from the positive antigen (the idea being that similar sequence
        TCRs bind similar antigens)
    """
    tmp_dfs = []

    pos_df_clean = pos_df.select(FORMAT_COLS + TCRDIST_COLS)
    neg_df_clean = negs_from_df.select(FORMAT_COLS + TCRDIST_COLS)

    pos_tcr_by_antigen = pos_df_clean.partition_by(
        FORMAT_ANTIGEN_COLS, maintain_order=True
    )

    neg_antigens = neg_df_clean.select(FORMAT_ANTIGEN_COLS).unique(maintain_order=True)

    for antigen_df in pos_tcr_by_antigen:

        # constant for every row in partition
        antigen = antigen_df.select(FORMAT_ANTIGEN_COLS)[0]

        pos_tcr = antigen_df.select(FORMAT_TCR_COLS + TCRDIST_COLS).unique(
            maintain_order=True
        )
        # in case focal antigen is in neg_antigens
        # below filtering covers this case, but just to be explicit
        other_antigens = neg_antigens.join(
            antigen, on=FORMAT_ANTIGEN_COLS, how="anti", maintain_order="left_right"
        )

        # filter out other_antigens based on pmhc distance

        peptide = antigen.select("peptide").item()
        mhc_1_seq = antigen.select("mhc_1_seq").item()
        mhc_2_seq = antigen.select("mhc_2_seq").item()

        other_antigens = (
            other_antigens.with_columns(
                pl.struct(
                    pl.col("peptide"),
                    pl.col("mhc_1_seq"),
                    pl.col("mhc_2_seq"),
                )
                .map_elements(
                    lambda x: (editdistance.eval(x["mhc_1_seq"], mhc_1_seq))
                    + (editdistance.eval(x["mhc_2_seq"], mhc_2_seq))
                    + (editdistance.eval(x["peptide"], peptide)),
                    return_dtype=pl.Int64,
                )
                .alias("edit_dist")
            )
            .filter(pl.col("edit_dist") > pmhc_edit_thresh)
            .drop("edit_dist")
        )

        # focal neg rows to sample TCRs from
        neg_rows = neg_df_clean.join(
            other_antigens, on=FORMAT_ANTIGEN_COLS, maintain_order="left_right"
        )
        neg_tcr = neg_rows.select(FORMAT_TCR_COLS + TCRDIST_COLS).unique(
            maintain_order=True
        )

        focal_rows = pl.concat([antigen_df, neg_rows])

        tcrdist_pos_neg, dist_mtx = pw_tcrdist(pl.concat([pos_tcr, neg_tcr]))

        focal_rows = focal_rows.join(
            tcrdist_pos_neg,
            on=FORMAT_TCR_COLS + TCRDIST_COLS,
            maintain_order="left_right",
        )

        other_antigen_tcr_indices_boolmask = np.zeros(dist_mtx.shape[0], dtype=np.bool)
        other_antigen_tcr_indices = np.sort(
            (
                focal_rows.join(
                    other_antigens,
                    on=FORMAT_ANTIGEN_COLS,
                    how="inner",
                    maintain_order="left_right",
                )
                .select("index")
                .unique(maintain_order=True)
            )
            .to_series()
            .to_numpy()
        )
        other_antigen_tcr_indices_boolmask[other_antigen_tcr_indices] = True

        focal_antigen_tcr_indices_boolmask = np.zeros(dist_mtx.shape[0], dtype=np.bool)
        focal_antigen_tcr_indices = np.sort(
            (
                focal_rows.join(
                    antigen,
                    on=FORMAT_ANTIGEN_COLS,
                    how="inner",
                    maintain_order="left_right",
                )
                .select("index")
                .unique(maintain_order=True)
            )
            .to_series()
            .to_numpy()
        )
        focal_antigen_tcr_indices_boolmask[focal_antigen_tcr_indices] = True
        # potentially valid noncognate TCRs: aren't within 120 tcrdist
        # of any of the focal antigens TCRs
        noncognate_tcr_indices = other_antigen_tcr_indices[
            np.all(
                dist_mtx[other_antigen_tcr_indices_boolmask][
                    :, focal_antigen_tcr_indices_boolmask
                ]
                > tcrdist_thresh,
                axis=1,
            )
        ]

        # a noncognate TCR can't be from the same antigen as a
        # TCR that is too close to a cognate TCR
        invalid_noncognate_tcr_indices = other_antigen_tcr_indices[
            np.any(
                dist_mtx[other_antigen_tcr_indices_boolmask][
                    :, focal_antigen_tcr_indices_boolmask
                ]
                <= tcrdist_thresh,
                axis=1,
            )
        ]

        tmp_idx_invalid = pl.DataFrame({"index": invalid_noncognate_tcr_indices})

        invalid_antigen = (
            focal_rows.join(tmp_idx_invalid, on="index", maintain_order="left_right")
            .select(FORMAT_ANTIGEN_COLS)
            .unique(maintain_order=True)
        )

        tmp_idx = pl.DataFrame({"index": noncognate_tcr_indices})

        # the second join here is redundant, but leaving it in case steps change
        # in the future
        tmp_df = (
            focal_rows.join(tmp_idx, on="index", maintain_order="left_right")
            .join(
                invalid_antigen,
                on=FORMAT_ANTIGEN_COLS,
                how="anti",
                maintain_order="left_right",
            )
            .select(FORMAT_TCR_COLS + TCRDIST_COLS + FORMAT_ANTIGEN_COLS)
            .unique(maintain_order=True)
        )

        # quick error check: no overlap in noncognate TCRs antigens with focal antigen
        if (
            tmp_df.join(
                antigen, on=FORMAT_ANTIGEN_COLS, maintain_order="left_right"
            ).height
            != 0
        ):
            raise ValueError

        # not suspect positive
        tmp_dfs.append(
            antigen.join(
                focal_rows.join(tmp_idx, on="index", maintain_order="left_right")
                .rename(SOURCE_RENAME_DICT)
                .select(FORMAT_TCR_COLS + TCRDIST_COLS + SOURCE_ANTIGEN_COLS)
                .group_by(FORMAT_TCR_COLS + SOURCE_ANTIGEN_COLS, maintain_order=True)
                .agg(pl.col(colname).first() for colname in TCRDIST_COLS)
                .unique(maintain_order=True),
                how="cross",
                maintain_order="left_right",
            )
        )

    noncognate = pl.concat(tmp_dfs, how="vertical")
    noncognate = generate_job_name(
        noncognate,
        [
            "peptide",
            "mhc_1_seq",
            "mhc_2_seq",
            "tcr_1_seq",
            "tcr_2_seq",
        ],
    )
    noncognate = noncognate.with_columns(pl.lit(False).alias("cognate"))
    noncognate = (
        noncognate.select(FORMAT_COLS + TCRDIST_COLS + SOURCE_ANTIGEN_COLS)
        .group_by(FORMAT_COLS, maintain_order=True)
        .agg(
            [pl.col(colname).first() for colname in TCRDIST_COLS] + SOURCE_ANTIGEN_COLS
        )
        .explode(SOURCE_ANTIGEN_COLS)
    ).select(FORMAT_COLS + TCRDIST_COLS + SOURCE_ANTIGEN_COLS)
    return noncognate


# def generate_all_possible_negs(
#     df,
#     tcrdist_thresh=120,
#     pmhc_edit_thresh=3,
#     cross_class=False,
# ):
#     """
#     DF must contain FORMAT_COLS as well as TCRDIST_COLS

#     Rows must be totally unqiue
#     """
#     tmp_dfs = []

#     all_antigens = df.select(FORMAT_ANTIGEN_COLS).unique(maintain_order=True)
#     all_tcrs = df.select(FORMAT_TCR_COLS + TCRDIST_COLS).unique(maintain_order=True)
#     all_tcrs, dist_mtx = pw_tcrdist(all_tcrs)

#     df_with_idx = df.join(
#         all_tcrs, on=FORMAT_TCR_COLS + TCRDIST_COLS, maintain_order="left_right"
#     )
#     found_set = set()

#     for row in (
#         df.select(FORMAT_ANTIGEN_COLS).unique(maintain_order=True).iter_rows(named=True)
#     ):

#         antigen = pl.DataFrame(row).select(FORMAT_ANTIGEN_COLS)

#         # tcrs_for_antigen = (
#         #     df.join(antigen, on=FORMAT_ANTIGEN_COLS, how="inner")
#         #     .select("index")
#         #     .unique(maintain_order=True)
#         #     .to_series()
#         #     .to_numpy()
#         # )

#         # don't favor TCRs from a particularly well-represented antigen
#         # do this by first randomly sampling an antigen,
#         # then choosing a random TCR from it

#         neg_found = 0

#         other_antigens = all_antigens.join(
#             antigen, on=FORMAT_ANTIGEN_COLS, how="anti", maintain_order="left_right"
#         )

#         if cross_class:
#             other_antigens = other_antigens.filter(
#                 pl.col("mhc_class")
#                 != antigen.select("mhc_class").unique(maintain_order=True)
#             )

#         # other_antigen = other_antigens.sample(n=1, shuffle=True)

#         # # do we reject based on antigen distance?

#         peptide = antigen.select("peptide").item()
#         mhc_1_seq = antigen.select("mhc_1_seq").item()
#         mhc_2_seq = antigen.select("mhc_2_seq").item()

#         other_antigens = other_antigens.with_columns(
#             pl.struct(
#                 pl.col("peptide"),
#                 pl.col("mhc_1_seq"),
#                 pl.col("mhc_2_seq"),
#             )
#             .map_elements(
#                 lambda x: (editdistance.eval(x["mhc_1_seq"], mhc_1_seq))
#                 + (editdistance.eval(x["mhc_2_seq"], mhc_2_seq))
#                 + (editdistance.eval(x["peptide"], peptide)),
#                 return_dtype=pl.Int64,
#             )
#             .alias("edit_dist")
#         ).filter(pl.col("edit_dist") > pmhc_edit_thresh)

#         # if antigen_distance <= pmhc_edit_thresh:
#         #     # reject, antigen too close
#         #     continue

#         # select a TCR
#         # tcrs_other_antigen = (
#         #     df_with_idx.join(
#         #         other_antigens, on=FORMAT_ANTIGEN_COLS, how="inner"
#         #     )
#         #     .select(FORMAT_TCR_COLS + TCRDIST_COLS + ["index"])
#         #     .unique(maintain_order=True)
#         # )

#         other_antigen_tcr_indices_boolmask = np.zeros(dist_mtx.shape[0], dtype=np.bool)
#         other_antigen_tcr_indices = np.sort(
#             (
#                 df_with_idx.join(
#                     other_antigens,
#                     on=FORMAT_ANTIGEN_COLS,
#                     how="inner",
#                     maintain_order="left_right",
#                 )
#                 .select("index")
#                 .unique(maintain_order=True)
#             )
#             .to_series()
#             .to_numpy()
#         )
#         other_antigen_tcr_indices_boolmask[other_antigen_tcr_indices] = True

#         # tcr_1_seq = tcr_other_antigen.select("tcr_1_seq").item()
#         # tcr_2_seq = tcr_other_antigen.select("tcr_2_seq").item()

#         # query = peptide + mhc_1_seq + mhc_2_seq + tcr_1_seq + tcr_2_seq

#         # # negs must be unqiue
#         # if query in found_set:
#         #     print("Miss")
#         #     continue

#         # target_row = tcr_other_antigen.select("index").item()

#         focal_antigen_tcr_indices_boolmask = np.zeros(dist_mtx.shape[0], dtype=np.bool)
#         focal_antigen_tcr_indices = (
#             (
#                 df_with_idx.join(
#                     antigen,
#                     on=FORMAT_ANTIGEN_COLS,
#                     how="inner",
#                     maintain_order="left_right",
#                 )
#                 .select("index")
#                 .unique(maintain_order=True)
#             )
#             .to_series()
#             .to_numpy()
#         )
#         focal_antigen_tcr_indices_boolmask[focal_antigen_tcr_indices] = True

#         noncognate_tcr_indices = other_antigen_tcr_indices[
#             np.all(
#                 dist_mtx[other_antigen_tcr_indices_boolmask][
#                     :, focal_antigen_tcr_indices_boolmask
#                 ]
#                 > tcrdist_thresh,
#                 axis=1,
#             )
#         ]

#         tmp_idx = pl.DataFrame({"index": noncognate_tcr_indices})

#         # # we only care about distance to binding TCRs
#         # # for the focal antigen
#         # tcr_hits = [
#         #     idx
#         #     for idx in np.argwhere(dist_mtx[target_row] <= tcrdist_thresh)[
#         #         :, 0
#         #     ]
#         #     if idx in focal_antigen_tcr_indices
#         # ]

#         # if len(tcr_hits) != 0:
#         #     print("Miss")
#         #     continue

#         tmp_df = (
#             df_with_idx.join(tmp_idx, on="index", maintain_order="left_right")
#             .select(FORMAT_TCR_COLS + TCRDIST_COLS + FORMAT_ANTIGEN_COLS)
#             .unique(maintain_order=True)
#         )

#         if (
#             tmp_df.join(
#                 antigen, on=FORMAT_ANTIGEN_COLS, maintain_order="left_right"
#             ).height
#             != 0
#         ):
#             raise ValueError

#         # not suspect positive
#         tmp_dfs.append(
#             antigen.join(
#                 df_with_idx.join(tmp_idx, on="index", maintain_order="left_right")
#                 .rename(SOURCE_RENAME_DICT)
#                 .select(FORMAT_TCR_COLS + TCRDIST_COLS + SOURCE_ANTIGEN_COLS)
#                 .group_by(FORMAT_TCR_COLS + SOURCE_ANTIGEN_COLS, maintain_order=True)
#                 .agg(pl.col(colname).first() for colname in TCRDIST_COLS)
#                 .unique(maintain_order=True),
#                 how="cross",
#                 maintain_order="left_right",
#             )
#         )
#         # neg_found += 1
#         # found_set.add(query)

#     noncognate = pl.concat(tmp_dfs, how="vertical")
#     noncognate = generate_job_name(
#         noncognate,
#         [
#             "peptide",
#             "mhc_1_seq",
#             "mhc_2_seq",
#             "tcr_1_seq",
#             "tcr_2_seq",
#         ],
#     )
#     noncognate = noncognate.with_columns(pl.lit(False).alias("cognate"))
#     noncognate = (
#         noncognate.select(FORMAT_COLS + TCRDIST_COLS + SOURCE_ANTIGEN_COLS)
#         .group_by(FORMAT_COLS, maintain_order=True)
#         .agg(
#             [pl.col(colname).first() for colname in TCRDIST_COLS] + SOURCE_ANTIGEN_COLS
#         )
#         .explode(SOURCE_ANTIGEN_COLS)
#     ).select(FORMAT_COLS + TCRDIST_COLS + SOURCE_ANTIGEN_COLS)
#     return noncognate


def generate_negatives_antigen_matched(
    df, tcrdist_thresh=120, pmhc_edit_thresh=3, n_neg=1, pos=None, dist=None
):
    """
    DF must contain FORMAT_COLS as well as TCRDIST_COLS

    Rows must be totally unqiue
    """

    tmp_dfs = []

    all_antigens = df.select(FORMAT_ANTIGEN_COLS).unique(maintain_order=True)

    if dist is None:
        all_tcrs = df.select(FORMAT_TCR_COLS + TCRDIST_COLS).unique(maintain_order=True)
        all_tcrs, dist_mtx = pw_tcrdist(all_tcrs)
        df_with_idx = df.join(
            all_tcrs, on=FORMAT_TCR_COLS + TCRDIST_COLS, maintain_order="left_right"
        )

    else:
        df_with_idx = df
        dist_mtx = dist

    found_set = set()

    iter_df = df.select(FORMAT_COLS) if pos is None else pos.select(FORMAT_COLS)

    for row in iter_df.iter_rows(named=True):

        antigen = pl.DataFrame(row).select(FORMAT_ANTIGEN_COLS)

        # tcrs_for_antigen = (
        #     df.join(antigen, on=FORMAT_ANTIGEN_COLS, how="inner")
        #     .select("index")
        #     .unique(maintain_order=True)
        #     .to_series()
        #     .to_numpy()
        # )

        # don't favor TCRs from a particularly well-represented antigen
        # do this by first randomly sampling an antigen,
        # then choosing a random TCR from it

        neg_found = 0

        while neg_found < n_neg:

            other_antigens = all_antigens.join(
                antigen, on=FORMAT_ANTIGEN_COLS, how="anti", maintain_order="left_right"
            )

            other_antigen = other_antigens.sample(n=1, shuffle=True)

            # do we reject based on antigen distance?

            peptide = antigen.select("peptide").item()
            mhc_1_seq = antigen.select("mhc_1_seq").item()
            mhc_2_seq = antigen.select("mhc_2_seq").item()

            antigen_distance = (
                other_antigen.with_columns(
                    pl.struct(
                        pl.col("peptide"),
                        pl.col("mhc_1_seq"),
                        pl.col("mhc_2_seq"),
                    )
                    .map_elements(
                        lambda x: (editdistance.eval(x["mhc_1_seq"], mhc_1_seq))
                        + (editdistance.eval(x["mhc_2_seq"], mhc_2_seq))
                        + (editdistance.eval(x["peptide"], peptide)),
                        return_dtype=pl.Int64,
                    )
                    .alias("edit_dist")
                )
                .select("edit_dist")
                .item()
            )

            if antigen_distance <= pmhc_edit_thresh:
                # reject, antigen too close
                continue

            # select a TCR
            tcr_other_antigen = (
                df_with_idx.join(
                    other_antigen,
                    on=FORMAT_ANTIGEN_COLS,
                    how="inner",
                    maintain_order="left_right",
                )
                .select(FORMAT_TCR_COLS + TCRDIST_COLS + ["index"])
                .unique(maintain_order=True)
            ).sample(n=1, shuffle=True)

            tcr_1_seq = tcr_other_antigen.select("tcr_1_seq").item()
            tcr_2_seq = tcr_other_antigen.select("tcr_2_seq").item()

            query = peptide + mhc_1_seq + mhc_2_seq + tcr_1_seq + tcr_2_seq

            # negs must be unqiue
            if query in found_set:
                # print("Miss")
                continue

            target_row = tcr_other_antigen.select("index").item()

            focal_antigen_tcr_indices = set(
                (
                    (
                        df_with_idx.join(
                            antigen,
                            on=FORMAT_ANTIGEN_COLS,
                            how="inner",
                            maintain_order="left_right",
                        )
                        .select("index")
                        .unique(maintain_order=True)
                    )
                    .to_series()
                    .to_list()
                )
            )

            # we only care about distance to binding TCRs
            # for the focal antigen
            tcr_hits = [
                idx
                for idx in np.argwhere(dist_mtx[target_row] <= tcrdist_thresh)[:, 0]
                if idx in focal_antigen_tcr_indices
            ]

            if len(tcr_hits) != 0:
                print("Miss")
                continue

            # not suspect positive
            tmp_dfs.append(
                antigen.join(
                    tcr_other_antigen.select(pl.exclude("index")),
                    how="cross",
                    maintain_order="left_right",
                )
            )
            neg_found += 1
            found_set.add(query)

    noncognate = pl.concat(tmp_dfs, how="vertical")
    noncognate = generate_job_name(
        noncognate,
        [
            "peptide",
            "mhc_1_seq",
            "mhc_2_seq",
            "tcr_1_seq",
            "tcr_2_seq",
        ],
    )
    noncognate = noncognate.with_columns(pl.lit(False).alias("cognate"))
    noncognate = noncognate.select(FORMAT_COLS + TCRDIST_COLS)
    return noncognate


def sample_to(antigen_df, neg_df):
    """
    For each antigen in antigen_df, select neg_per_pos negatives from neg_df.

    antigen_df must have a column "needed_negs"

    Negative selection follows these rules:

    1. If there aren't enough negatives in neg_df, add all of them, and add
        the antigen to the returned missing_df with a needed_negatives column
    2. If there are enough negatives in neg_df, select neg_per_pos negatives.
        Select using the sampling method:
        1. Select a random source antigen
        2. Select a random TCR from that source antigen
        3. Remove that TCR from all source antigens (since they may be shared across source antigens)
    """

    new_neg_df = []
    missing_df = []

    for row in tqdm(
        antigen_df.iter_rows(named=True),
        total=antigen_df.height,
        desc="Processing rows",
    ):
        antigen = pl.DataFrame([row])
        desired_neg = row["needed_negs"]

        # remove source, which is duplicated
        true_neg_df = (
            neg_df.join(antigen, on=FORMAT_ANTIGEN_COLS, maintain_order="left_right")
            .select(FORMAT_COLS + TCRDIST_COLS)
            .unique(maintain_order=True)
        )
        true_neg = true_neg_df.height

        if desired_neg > true_neg:
            needed_neg = desired_neg - true_neg
            # jsut add them all
            new_neg_df.append(true_neg_df)
            missing_df.append(
                antigen.with_columns(pl.lit(needed_neg).alias("needed_negs"))
            )

        elif true_neg > desired_neg:

            candidate_neg = neg_df.join(
                antigen, on=FORMAT_ANTIGEN_COLS, maintain_order="left_right"
            )
            neg_by_source_antigen = candidate_neg.partition_by(
                SOURCE_ANTIGEN_COLS, maintain_order=True
            )
            neg_partition_meta = [df.height for df in neg_by_source_antigen]
            remaining_neg = list(range(len(neg_by_source_antigen)))

            selected = 0
            sel_tcrs = []

            while selected < desired_neg:

                ant_int = random.choice(remaining_neg)

                # neg_partition_meta[ant_int] -= 1
                selected += 1

                # sample a TCR and strike it out from all source antigen
                sel_tcr = (
                    neg_by_source_antigen[ant_int]
                    .select(FORMAT_TCR_COLS + TCRDIST_COLS)
                    .unique(maintain_order=True)
                    .sample(n=1, shuffle=True)
                )
                sel_tcrs.append(sel_tcr)

                rmvlist = []
                for i in remaining_neg:
                    # strikeout TCR
                    neg_by_source_antigen[i] = neg_by_source_antigen[i].join(
                        sel_tcr,
                        on=FORMAT_TCR_COLS + TCRDIST_COLS,
                        how="anti",
                        maintain_order="left_right",
                    )

                    # update remaining TCRs to select
                    neg_partition_meta[i] = neg_by_source_antigen[i].height

                    if neg_partition_meta[i] == 0:
                        rmvlist.append(i)

                for i in rmvlist:
                    remaining_neg.remove(i)

            sel_tcrs = pl.concat(sel_tcrs)

            noncognate = antigen.join(
                sel_tcrs, how="cross", maintain_order="left_right"
            )
            noncognate = generate_job_name(
                noncognate,
                [
                    "peptide",
                    "mhc_1_seq",
                    "mhc_2_seq",
                    "tcr_1_seq",
                    "tcr_2_seq",
                ],
            )
            noncognate = noncognate.with_columns(pl.lit(False).alias("cognate")).select(
                FORMAT_COLS + TCRDIST_COLS
            )

            if noncognate.height != noncognate.unique(maintain_order=True).height:
                print(row)

            new_neg_df.append(noncognate)

        else:
            # exactly the right number of negs
            new_neg_df.append(true_neg_df)

    if len(missing_df) == 0:
        missing = None
    else:
        missing = pl.concat(missing_df, how="vertical_relaxed")

    return pl.concat(new_neg_df, how="vertical_relaxed"), missing


# def sample_supplemental_negatives(missing_neg_antigens, cross_class_negatives):
#     """
#     Takes the missing_neg_antigens df from sample_to and the cross_class_negatives
#     which is a dataframe of all possible negative TCRs from the opposite MHC class.

#     For each antigen in missing_neg_antigens, select the number of negatives
#     needed from the cross_class_negatives. Use the same sampling method as
#     sample_to, but only select from the cross_class_negatives.
#     """
#     new_neg_df = []
#     missing_df = []

#     for row in tqdm(
#         missing_neg_antigens.iter_rows(named=True),
#         total=missing_neg_antigens.height,
#         desc="Processing rows",
#     ):
#         antigen = pl.DataFrame([row])

#         needed_neg = row["needed_negatives"]

#         candidate_neg = cross_class_negatives.join(
#             antigen, on=FORMAT_ANTIGEN_COLS, maintain_order="left_right"
#         )
#         neg_by_source_antigen = candidate_neg.partition_by(
#             SOURCE_ANTIGEN_COLS, maintain_order=True
#         )
#         neg_partition_meta = [df.height for df in neg_by_source_antigen]
#         remaining_neg = list(range(len(neg_by_source_antigen)))

#         if needed_neg > candidate_neg.height:
#             raise ValueError("Not enough negatives")

#         selected = 0
#         sel_tcrs = []

#         while selected < needed_neg:

#             ant_int = random.choice(remaining_neg)

#             # neg_partition_meta[ant_int] -= 1
#             selected += 1

#             # sample a TCR and strike it out from all source antigen
#             sel_tcr = (
#                 neg_by_source_antigen[ant_int]
#                 .select(FORMAT_TCR_COLS + TCRDIST_COLS)
#                 .unique(maintain_order=True)
#                 .sample(n=1, shuffle=True)
#             )
#             sel_tcrs.append(sel_tcr)

#             rmvlist = []
#             for i in remaining_neg:
#                 # strikeout TCR
#                 neg_by_source_antigen[i] = neg_by_source_antigen[i].join(
#                     sel_tcr,
#                     on=FORMAT_TCR_COLS + TCRDIST_COLS,
#                     how="anti",
#                     maintain_order="left_right",
#                 )

#                 # update remaining TCRs to select
#                 neg_partition_meta[i] = neg_by_source_antigen[i].height

#                 if neg_partition_meta[i] == 0:
#                     rmvlist.append(i)

#             for i in rmvlist:
#                 remaining_neg.remove(i)

#         sel_tcrs = pl.concat(sel_tcrs)

#         noncognate = antigen.join(sel_tcrs, how="cross", maintain_order="left_right")
#         noncognate = generate_job_name(
#             noncognate,
#             [
#                 "peptide",
#                 "mhc_1_seq",
#                 "mhc_2_seq",
#                 "tcr_1_seq",
#                 "tcr_2_seq",
#             ],
#         )
#         noncognate = noncognate.with_columns(pl.lit(False).alias("cognate")).select(
#             FORMAT_COLS + TCRDIST_COLS
#         )

#         if noncognate.height != noncognate.unique(maintain_order=True).height:
#             print(row)

#         new_neg_df.append(noncognate)

#     return pl.concat(new_neg_df, how="vertical_relaxed")
