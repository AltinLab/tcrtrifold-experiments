
from mdaf3.AF3OutputParser import AF3Output
from mdaf3.FeatureExtraction import *
from tcr_format_parsers.common.TCRUtils import (
    annotate_tcr,
    tcr_by_imgt_region,
    pw_tcrdist,
)
from tcr_format_parsers.common.TriadUtils import (
    FORMAT_TCR_COLS,
    FORMAT_ANTIGEN_COLS,
    FORMAT_MHC_COLS,
    TCRDIST_COLS,
)
import numpy as np
import polars as pl
from scipy.stats import gmean
import scipy
from itertools import repeat


TCRDOCK_COLS = [
    "d",
    "torsion",
    "tcr_unit_y",
    "tcr_unit_z",
    "mhc_unit_y",
    "mhc_unit_z",
]


def tcrdock_info_as_feat(row, dockinfo_path, fname_col="job_name"):
    if not (dockinfo_path / (row[fname_col] + ".tsv")).exists():
        for key in TCRDOCK_COLS:
            row[key] = None
        return row

    pred_dockinfo = pl.read_csv(
        dockinfo_path / (row[fname_col] + ".tsv"), separator="\t"
    )

    pred_dockinfo = pred_dockinfo.to_dicts()[0]

    for key in TCRDOCK_COLS:
        row[key] = pred_dockinfo[key]

    return row


def extract_mean_tcr_pmhc_pae(row, af3_parent_dir, **kwargs):
    af3_output = AF3Output(af3_parent_dir / row["job_name"], **kwargs)
    u = af3_output.get_mda_universe(**kwargs)

    peptide_res = u.select_atoms("segid A").residues

    if row["mhc_class"] == "II":
        mhc_residx = u.select_atoms("segid B or segid C").residues.resindices
    else:
        mhc_residx = u.select_atoms("segid B").residues.resindices

    tcr_residx = u.select_atoms("segid D or segid E").residues.resindices

    pae = af3_output.get_pae_ndarr(**kwargs)

    if row["mhc_class"] == "II":
        if len(peptide_res) <= 9:
            row["mean_p_tcr_pae"] = (
                np.mean(pae[peptide_res.resindices][:, tcr_residx]) / 100
            )
            row["mean_tcr_p_pae"] = (
                np.mean(pae[tcr_residx][:, peptide_res.resindices]) / 100
            )
        else:

            p_tcr_window_means = []
            tcr_p_window_means = []
            for i in range(len(peptide_res) - 8):
                p_tcr_window_means.append(
                    (
                        np.mean(
                            pae[peptide_res[i : i + 9].resindices][
                                :, tcr_residx
                            ]
                        )
                        / 100
                    )
                )
                tcr_p_window_means.append(
                    (
                        np.mean(
                            pae[tcr_residx][
                                :, peptide_res[i : i + 9].resindices
                            ]
                        )
                        / 100
                    )
                )
            row["mean_p_tcr_pae"] = gmean(p_tcr_window_means)
            row["mean_tcr_p_pae"] = gmean(tcr_p_window_means)

    else:
        row["mean_p_tcr_pae"] = (
            np.mean(pae[peptide_res.resindices][:, tcr_residx]) / 100
        )
        row["mean_tcr_p_pae"] = (
            np.mean(pae[tcr_residx][:, peptide_res.resindices]) / 100
        )

    # mean_interface_pae = np.mean(pae[pmhc_residx][:, tcr_residx])

    row["mean_mhc_tcr_pae"] = np.mean(pae[mhc_residx][:, tcr_residx])
    row["mean_tcr_mhc_pae"] = np.mean(pae[tcr_residx][:, mhc_residx])

    return row


def extract_mean_peptide_mhc_pae(row, af3_parent_dir, **kwargs):
    af3_output = AF3Output(af3_parent_dir / row["job_name"], **kwargs)
    u = af3_output.get_mda_universe(**kwargs)

    peptide_res = u.select_atoms("segid A").residues

    if row["mhc_class"] == "II":
        mhc_residx = u.select_atoms("segid B or segid C").residues.resindices
    else:
        mhc_residx = u.select_atoms("segid B").residues.resindices

    pae = af3_output.get_pae_ndarr(**kwargs)

    if row["mhc_class"] == "II":
        if len(peptide_res) <= 9:
            row["mean_p_mhc_pae"] = (
                np.mean(pae[peptide_res.resindices][:, mhc_residx]) / 100
            )
        else:

            p_mhc_window_means = []
            for i in range(len(peptide_res) - 8):
                p_mhc_window_means.append(
                    (
                        np.mean(
                            pae[peptide_res[i : i + 9].resindices][
                                :, mhc_residx
                            ]
                        )
                        / 100
                    )
                )

            row["mean_p_mhc_pae"] = gmean(p_mhc_window_means)

    else:
        row["mean_p_mhc_pae"] = (
            np.mean(pae[peptide_res.resindices][:, mhc_residx]) / 100
        )

    return row


def extract_min_tcr_pmhc_pae(row, af3_parent_dir, **kwargs):
    af3_output = AF3Output(af3_parent_dir / row["job_name"], **kwargs)
    summ = af3_output.get_summary_metrics(**kwargs)

    pae_chain = np.array(summ["chain_pair_pae_min"])
    p_tcr = min(pae_chain[0][3], pae_chain[0][4])
    tcr_p = min(pae_chain[3][0], pae_chain[4][0])

    if row["mhc_class"] == "II":
        # min of both mhc chains
        m_tcr = min(
            min(pae_chain[1][3], pae_chain[1][4]),
            min(pae_chain[2][3], pae_chain[2][4]),
        )
        tcr_m = min(
            min(pae_chain[3][1], pae_chain[4][1]),
            min(pae_chain[3][2], pae_chain[4][2]),
        )
    else:
        # min of chain 1
        m_tcr = min(
            pae_chain[1][3],
            pae_chain[1][4],
        )
        tcr_m = min(
            pae_chain[3][1],
            pae_chain[4][1],
        )
    min_pae_dict = {
        "min_p_tcr_pae": p_tcr,
        "min_mhc_tcr_pae": m_tcr,
        "min_tcr_p_pae": tcr_p,
        "min_tcr_mhc_pae": tcr_m,
    }

    row.update(min_pae_dict)

    return row


def extract_summary_metrics(row, af3_parent_dir, **kwargs):
    af3_output = AF3Output(af3_parent_dir / row["job_name"], **kwargs)

    summ = af3_output.get_summary_metrics(**kwargs)

    # scalar_dict = {
    #     "ptm": summ["ptm"],
    #     "iptm": summ["iptm"],
    #     "fraction_disordered": summ["fraction_disordered"],
    #     "has_clash": summ["has_clash"],
    #     "ranking_score": summ["ranking_score"],
    # }

    small_arrs = [
        "chain_pair_pae_min",
        "chain_pair_iptm",
        "chain_ptm",
        "chain_iptm",
    ]

    # convert to list of lists
    for arr in small_arrs:
        summ[arr] = summ[arr].tolist()

    row.update(summ)

    # schema_overrides = {
    #     "chain_pair_pae_min": pl.List(pl.List(pl.Float32)),
    #     "chain_pair_iptm": pl.List(pl.List(pl.Float32)),
    #     "chain_ptm": pl.List(pl.Float32),
    #     "chain_iptm": pl.List(pl.Float32),
    # }

    return row


def extract_peptide_pLDDT(row, af3_parent_dir):
    af3_output = AF3Output(af3_parent_dir / row["job_name"])
    u_af3 = af3_output.get_mda_universe()

    peptide_sel = u_af3.select_atoms("segid A").residues

    if row["mhc_class"] == "II":

        if len(peptide_sel) <= 9:
            row["peptide_mean_pLDDT"] = 1 - (
                peptide_sel.atoms.tempfactors.mean() / 100
            )

        else:

            window_means = []
            for i in range(len(peptide_sel) - 8):
                window_means.append(
                    1 - (peptide_sel[i : i + 9].atoms.tempfactors.mean() / 100)
                )

            row["peptide_mean_pLDDT"] = gmean(window_means)

    else:
        row["peptide_mean_pLDDT"] = 1 - (
            peptide_sel.atoms.tempfactors.mean() / 100
        )

    return row


def extract_cdr_pLDDT(row, af3_parent_dir):
    af3_output = AF3Output(af3_parent_dir / row["job_name"])
    u_af3 = af3_output.get_mda_universe()

    for segid, tcr_num in zip(["D", "E"], [1, 2]):

        tcr_resindices = u_af3.select_atoms(
            f"segid {segid}"
        ).residues.resindices

        tcr_indices, imgt_num, _ = annotate_tcr(
            row[f"tcr_{tcr_num}_seq"],
            tcr_resindices,
            row[f"tcr_{tcr_num}_chain"],
            row[f"tcr_{tcr_num}_species"],
        )

        cdr_1_pLDDT = u_af3.residues[
            tcr_indices[((imgt_num >= 27) & (imgt_num <= 38))]
        ].atoms.tempfactors

        row[f"tcr_{tcr_num}_cdr_1_mean_pLDDT"] = cdr_1_pLDDT.mean()

        cdr_2_pLDDT = u_af3.residues[
            tcr_indices[((imgt_num >= 56) & (imgt_num <= 65))]
        ].atoms.tempfactors

        row[f"tcr_{tcr_num}_cdr_2_mean_pLDDT"] = cdr_2_pLDDT.mean()

        cdr_2_5_pLDDT = u_af3.residues[
            tcr_indices[((imgt_num >= 81) & (imgt_num <= 86))]
        ].atoms.tempfactors

        row[f"tcr_{tcr_num}_cdr_2_5_mean_pLDDT"] = cdr_2_5_pLDDT.mean()

        cdr_3_pLDDT = u_af3.residues[
            tcr_indices[((imgt_num >= 104) & (imgt_num <= 118))]
        ].atoms.tempfactors

        row[f"tcr_{tcr_num}_cdr_3_mean_pLDDT"] = cdr_3_pLDDT.mean()

    return row


def extract_mhc_helix_pLDDT(row, af3_parent_dir):
    af3_output = AF3Output(af3_parent_dir / row["job_name"])
    u_af3 = af3_output.get_mda_universe()

    if row["mhc_class"] == "I":
        mhc_atoms = u_af3.select_atoms("segid B")
    else:
        mhc_atoms = u_af3.select_atoms("segid B or segid C")

    raw_helix_ix = raw_helix_indices(mhc_atoms)

    helix_resindices = join_helix_indices(raw_helix_ix)

    if len(helix_resindices) != 2:
        raise ValueError(
            f"Expected 2 helices for MHC-{row['mhc_class']} chain, "
            f"got {len(helix_resindices)}"
        )

    # flatten list of lists
    helix_resindices = helix_resindices[0] + helix_resindices[1]

    helix_pLDDT = u_af3.residues[helix_resindices].atoms.tempfactors

    row["mhc_helices_mean_pLDDT"] = helix_pLDDT.mean()

    return row


def extract_num_contacts(row, af3_parent_dir, **kwargs):
    af3_output = AF3Output(af3_parent_dir / row["job_name"], **kwargs)
    u_af3 = af3_output.get_mda_universe()
    contact_probs = af3_output.get_contact_prob_ndarr()

    peptide_tcr_thresh = 0.5
    tcr_mhc_thresh = 0.9

    if row["mhc_class"] == "I":
        mhc_atoms = u_af3.select_atoms("segid B")

        raw_helix_ix = raw_helix_indices(mhc_atoms)
        helix_resindices = join_helix_indices(raw_helix_ix)

        hla_a_res = helix_resindices[0]
        hla_b_res = helix_resindices[1]

    else:
        hla_a_res = u_af3.select_atoms("segid B").residues.resindices
        hla_b_res = u_af3.select_atoms("segid C").residues.resindices

    # first, extract broad summary contact metrics
    peptide_sel = u_af3.select_atoms("segid A").residues
    peptide_res_all = peptide_sel.resindices
    tcr_res = u_af3.select_atoms("segid D or segid E").residues.resindices

    row["tcr_mhc_contacts"] = np.count_nonzero(
        contact_probs[tcr_res][:, hla_a_res] > tcr_mhc_thresh
    ) + np.count_nonzero(contact_probs[tcr_res][:, hla_b_res] > tcr_mhc_thresh)

    row["peptide_tcr_contacts"] = np.count_nonzero(
        contact_probs[tcr_res][:, peptide_res_all] > peptide_tcr_thresh
    )

    # find best 9-mer
    max_min = 0
    pep_start = None

    if len(peptide_sel) < 9:
        row["contact_map"] = None
        return row

    for i in range(len(peptide_sel) - 8):
        min_pLDDT = peptide_sel[i : i + 9].atoms.tempfactors.min()

        if min_pLDDT > max_min:
            max_min = min_pLDDT
            pep_start = i

    pep_resindices = peptide_sel[pep_start : pep_start + 9].resindices

    tcr_1_resdict = tcr_by_imgt_region(
        row["tcr_1_seq"],
        u_af3.select_atoms("segid D").residues.resindices,
        row["tcr_1_chain"],
        row["tcr_1_species"],
    )

    tcr_2_resdict = tcr_by_imgt_region(
        row["tcr_2_seq"],
        u_af3.select_atoms("segid E").residues.resindices,
        row["tcr_2_chain"],
        row["tcr_2_species"],
    )

    contact_map = np.zeros((14, 11), dtype=np.float32)

    keys = ["fwr_1", "cdr_1", "fwr_2", "cdr_2", "fwr_3", "cdr_3", "fwr_4"]

    for i, key in enumerate(keys):

        tcr_1_rowidx = i
        tcr_2_rowidx = i + len(keys)

        for j in range(len(pep_resindices)):

            contact_map[tcr_1_rowidx][j] = np.count_nonzero(
                contact_probs[tcr_1_resdict[key]][:, pep_resindices[j]]
                > peptide_tcr_thresh
            )

            contact_map[tcr_2_rowidx][j] = np.count_nonzero(
                contact_probs[tcr_2_resdict[key]][:, pep_resindices[j]]
                > peptide_tcr_thresh
            )

        contact_map[tcr_1_rowidx][9] = np.count_nonzero(
            contact_probs[tcr_1_resdict[key]][:, hla_a_res] > tcr_mhc_thresh
        )

        contact_map[tcr_2_rowidx][9] = np.count_nonzero(
            contact_probs[tcr_2_resdict[key]][:, hla_a_res] > tcr_mhc_thresh
        )

        contact_map[tcr_1_rowidx][10] = np.count_nonzero(
            contact_probs[tcr_1_resdict[key]][:, hla_b_res] > tcr_mhc_thresh
        )

        contact_map[tcr_2_rowidx][10] = np.count_nonzero(
            contact_probs[tcr_2_resdict[key]][:, hla_b_res] > tcr_mhc_thresh
        )

    row["contact_map"] = contact_map.tolist()

    return row


def extract_imgt_correct_species(row, af3_parent_dir):

    for tcr_num in [1, 2]:
        try:
            annotate_tcr(
                row[f"tcr_{tcr_num}_seq"],
                # dummy arg
                np.arange(len(row[f"tcr_{tcr_num}_seq"])),
                row[f"tcr_{tcr_num}_chain"],
                row[f"tcr_{tcr_num}_species"],
                strict=True,
            )
        except ValueError as e:
            row[f"tcr_{tcr_num}_species_correct"] = False

        else:
            row[f"tcr_{tcr_num}_species_correct"] = True

    return row


def extract_tcr_pmhc_iptm_perseed(row, af3_parent_dir):
    af3_output = AF3Output(af3_parent_dir / row["job_name"])

    with af3_output._get_h5_handle() as hf:
        ranking_score_np = hf["ranking_scores"]["ranking_score"][:]
        seed_np = hf["ranking_scores"]["seed"][:]

        ranking_argsort_desc = np.argsort(ranking_score_np)[::-1]

        seeds_ranked = seed_np[ranking_argsort_desc].tolist()

    seeds_rank_dict = {num: seeds_ranked.index(num) for num in range(1, 6)}

    small_arrs = [
        "chain_pair_pae_min",
        "chain_pair_iptm",
    ]
    mhc_tcr_iptm_feats = []
    p_tcr_iptm_feats = []
    mhc_tcr_pae_feats = []
    p_tcr_pae_feats = []

    seeds_ran = []
    # convert to list of lists
    for seed in [1, 2, 3, 4, 5]:
        seeds_ran.append(seed)
        pl_index = 6
        best_seed = None
        for seed_ran in seeds_ran:
            if seeds_rank_dict[seed_ran] < pl_index:
                best_seed = seed_ran
                pl_index = seeds_rank_dict[seed_ran]

        # now extract features for curr best seed
        summ = af3_output.get_summary_metrics(seed=best_seed)
        pae = summ["chain_pair_pae_min"]
        iptm = summ["chain_pair_iptm"]
        row["mhc_tcr_iptm_" + str(seed)] = np.mean(iptm[1:3, 3:5])
        row["p_tcr_iptm_" + str(seed)] = np.mean(iptm[0, 3:5])
        row["mhc_tcr_pae_" + str(seed)] = np.mean(pae[1:3, 3:5])
        row["p_tcr_pae_" + str(seed)] = np.mean(pae[0, 3:5])

    return row


def raw_helix_indices(sel):
    # find helices
    # https://docs.mdanalysis.org/2.8.0/documentation_pages/analysis/dssp.html
    helix_resindices_boolmask = DSSP(sel).run().results.dssp_ndarray[0, :, 1]
    return sel.residues[helix_resindices_boolmask].resindices


def join_helix_indices(helix_resindices, len_tol=15, gap_tol=10):

    helices = []
    curr_helix = []
    prev_ix = helix_resindices[0] - 1
    for i in range(len(helix_resindices)):

        if helix_resindices[i] - prev_ix <= gap_tol:
            # fill gap
            append_ix = range(prev_ix + 1, helix_resindices[i] + 1)
            curr_helix.extend(append_ix)

            if i == len(helix_resindices) - 1 and len(curr_helix) >= len_tol:
                helices.append(curr_helix)

        else:
            if len(curr_helix) >= len_tol:
                helices.append(curr_helix)
            # start new helix
            curr_helix = [helix_resindices[i]]

        prev_ix = helix_resindices[i]

    return helices


def transform_zero_one(df, featnames_dict):

    featnames = list(featnames_dict.keys())
    for feat, direct_relationship in featnames_dict.items():

        max_feat = df.select(pl.col(feat).max()).item()
        min_feat = df.select(pl.col(feat).min()).item()

        if direct_relationship == True:
            df = df.with_columns(
                ((pl.col(feat) - min_feat) / (max_feat - min_feat)).alias(feat)
            )

        elif direct_relationship == False:

            df = df.with_columns(
                (
                    1 - ((pl.col(feat) - min_feat) / (max_feat - min_feat))
                ).alias(feat)
            )

    return df


def effective_variance(df):

    tmp_dfs = []
    for antigen in (
        df.select(FORMAT_ANTIGEN_COLS).unique().iter_rows(named=True)
    ):
        focal_tcr_cognate = df.join(
            pl.DataFrame(antigen), on=FORMAT_ANTIGEN_COLS
        ).filter(pl.col("cognate"))
        feats = cossin_embed(focal_tcr_cognate.select(TCRDOCK_COLS).to_numpy())
        R = np.corrcoef(feats, rowvar=False)

        eigs = np.linalg.eigvalsh(R)
        gen_var = np.prod(eigs)

        antigen_w_var = antigen.copy()
        antigen_w_var["tot_var_cognate"] = np.sum(np.cov(feats, rowvar=False))

        focal_tcr_noncognate = df.join(
            pl.DataFrame(antigen), on=FORMAT_ANTIGEN_COLS
        ).filter(~pl.col("cognate"))
        feats = cossin_embed(
            focal_tcr_noncognate.select(TCRDOCK_COLS).to_numpy()
        )
        antigen_w_var["tot_var_noncognate"] = np.sum(
            np.cov(feats, rowvar=False)
        )

        tmp_dfs.append(pl.DataFrame(antigen_w_var))

    return pl.concat(tmp_dfs)


def geomean_combine_confidence(
    df, featnames_dict, featnames_ranges, list_cols=[]
):
    # carefully transform each feature so that its
    # best (highest P(cognate)) value is 0 and its worst
    # value is 1

    # we will use the actual max and min ranges rather than
    # observed ranges in the dataset
    df_inv = df

    for featname, rel in featnames_dict.items():
        ra = featnames_ranges[featname]

        # direct relationship, higher = better prediction
        if rel:
            # invert
            df_inv = df_inv.with_columns(
                (1 - (pl.col(featname) / ra[1])).alias(featname)
            )

        else:
            df_inv = df_inv.with_columns(
                (pl.col(featname) / ra[1]).alias(featname)
            )

    df_agg = df_inv.group_by("entity_id").agg(
        [
            pl.col(colname).first().alias(colname)
            for colname in FORMAT_MHC_COLS
            + FORMAT_TCR_COLS
            + TCRDIST_COLS
            + ["group"]
            + ["cognate"]
        ]
        + [
            pl.col(colname).drop_nulls().flatten().unique()
            for colname in ["references", "assay_type", "receptor_id"]
        ]
        + [pl.col(colname) for colname in list_cols]
        + [
            # geometric mean
            # then transform so higher = better
            (1 - (pl.col(colname).log().mean().exp())).alias(colname)
            for colname in featnames_dict.keys()
        ]
    )

    return df_agg


def per_antigen_tcrdist_clust(df, use_provided_cdr=False):

    df_by_antigen = df.partition_by(
        FORMAT_ANTIGEN_COLS,
    )

    out_dfs = []

    for antigen_df in df_by_antigen:

        cdr_2_5_col = (
            ["tcr_1_cdr_2_5", "tcr_2_cdr_2_5"] if use_provided_cdr else []
        )
        tcr_df = antigen_df.select(
            FORMAT_TCR_COLS + TCRDIST_COLS + cdr_2_5_col
        )

        tcr_with_idx, pw_dist = pw_tcrdist(
            tcr_df,
            use_provided_cdr=use_provided_cdr,
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
        tcr_with_idx = tcr_with_idx.with_columns(
            pl.Series("cluster", clusters)
        )

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


def pw_tcrgeom_dist(tcr_df):

    tcr_geom_np_1 = cossin_embed(tcr_df.select(TCRDOCK_COLS).to_numpy())
    tcr_geom_np_2 = tcr_geom_np_1.copy()

    cognate_geom = cossin_embed(
        tcr_df.filter(pl.col("cognate")).select(TCRDOCK_COLS).to_numpy()
    )

    mu = cognate_geom.mean(axis=0)
    cov = np.cov(cognate_geom, rowvar=False)
    inv_cov = np.linalg.inv(cov)

    pw_dist = scipy.spatial.distance.cdist(
        tcr_geom_np_1, tcr_geom_np_2, metric="mahalanobis", VI=inv_cov
    )

    return pw_dist


def tcrdist_tcrdock_geomdist_corr(df, cognate_only=False):

    df_by_antigen = df.partition_by(
        FORMAT_ANTIGEN_COLS,
    )

    out_dfs = []

    for antigen_df in df_by_antigen:

        if cognate_only:
            # only keep cognate TCRs
            antigen_df = antigen_df.filter(pl.col("cognate"))

        tcr_df = antigen_df.select(
            ["cognate"] + FORMAT_TCR_COLS + TCRDOCK_COLS + TCRDIST_COLS
        )

        tcr_with_idx, pw_tcrdist_np = pw_tcrdist(
            tcr_df,
        )

        pw_tcrdock_geom_dist_np = pw_tcrgeom_dist(tcr_with_idx)

        compressed_geom_dist = scipy.spatial.distance.squareform(
            pw_tcrdock_geom_dist_np
        )
        compressed_seq_dist = scipy.spatial.distance.squareform(pw_tcrdist_np)

        r, p = spearmanr(compressed_geom_dist, compressed_seq_dist)

        antigen_df_corr = (
            antigen_df.select(FORMAT_ANTIGEN_COLS)
            .unique()
            .with_columns(
                [
                    pl.lit(r).alias("r"),
                    pl.lit(p).alias("p"),
                ]
            )
        )

        out_dfs.append(antigen_df_corr)

    return pl.concat(out_dfs)


from MDAnalysis.analysis import align, rms


def _pw_pep_rmsd(j1, j2, inf_path):
    af3_o_1 = AF3Output(inf_path / j1)
    af3_o_2 = AF3Output(inf_path / j2)

    u_1 = af3_o_1.get_mda_universe()
    u_2 = af3_o_2.get_mda_universe()

    align.alignto(u_1, u_2, select="segid B")

    return rms.rmsd(
        u_1.select_atoms("segid A and name CA").positions,
        u_2.select_atoms("segid A and name CA").positions,
        center=False,
        superposition=False,
    )


def pw_pep_rmsd(tcr_df, inf_path):

    job_names = tcr_df.select("job_name").to_series().to_list()

    pw_matrix = np.array([[(j1, j2) for j2 in job_names] for j1 in job_names])
    pw_condensed = pw_matrix[np.triu_indices(len(job_names), k=1)].tolist()
    rmsd_condensed = []

    # for job_1, job_2 in pw_condensed:
    #     rmsd_condensed.append(_pw_pep_rmsd(job_1, job_2, inf_path))
    jobs1, jobs2 = zip(*pw_condensed)
    rmsd_condensed = process_map(
        _pw_pep_rmsd,
        jobs1,
        jobs2,
        repeat(inf_path),
        chunksize=15,
        total=len(jobs1),
    )

    return np.array(rmsd_condensed)


def tcrdist_pep_rmsd_corr(df, inf_path, cognate_only=False):

    df_by_antigen = sorted(
        df.partition_by(
            FORMAT_ANTIGEN_COLS,
        ),
        key=lambda x: x.height,
    )

    out_dfs = []

    for antigen_df in df_by_antigen:

        if cognate_only:
            # only keep cognate TCRs
            antigen_df = antigen_df.filter(pl.col("cognate"))

        tcr_df = antigen_df.select(
            ["job_name", "cognate"]
            + FORMAT_TCR_COLS
            + TCRDIST_COLS
            + FORMAT_ANTIGEN_COLS
        )

        tcr_with_idx, pw_tcrdist_np = pw_tcrdist(
            tcr_df,
        )

        compressed_rmsd = pw_pep_rmsd(tcr_with_idx, inf_path)
        compressed_seq_dist = scipy.spatial.distance.squareform(pw_tcrdist_np)

        r, p = spearmanr(compressed_rmsd, compressed_seq_dist)

        print(f"pep {antigen_df.select("peptide")[0].item()}: r: {r}, p: {p}")

        antigen_df_corr = (
            antigen_df.select(FORMAT_ANTIGEN_COLS)
            .unique()
            .with_columns(
                [
                    pl.lit(r).alias("r"),
                    pl.lit(p).alias("p"),
                ]
            )
        )

        out_dfs.append(antigen_df_corr)

    return pl.concat(out_dfs)


extract_funcs = [
    extract_mean_tcr_pmhc_pae,
    extract_min_tcr_pmhc_pae,
    extract_summary_metrics,
    extract_peptide_pLDDT,
    extract_cdr_pLDDT,
    extract_mhc_helix_pLDDT,
    extract_num_contacts,
    extract_imgt_correct_species,
    extract_tcr_pmhc_iptm_perseed,
]

antigen_extract_funcs = [
    extract_mean_peptide_mhc_pae,
    extract_summary_metrics,
    extract_peptide_pLDDT,
    extract_mhc_helix_pLDDT,
]

all_featnames_dict = {
    "mean_tcr_p_pae": False,
    "mean_p_tcr_pae": False,
    "mean_mhc_tcr_pae": False,
    "mean_tcr_mhc_pae": False,
    "min_p_tcr_pae": False,
    "min_mhc_tcr_pae": False,
    "min_tcr_p_pae": False,
    "min_tcr_mhc_pae": False,
    "ptm": True,
    "iptm": True,
    "fraction_disordered": True,
    "ranking_score": True,
    "peptide_mean_pLDDT": False,
    "tcr_1_cdr_1_mean_pLDDT": True,
    "tcr_1_cdr_2_mean_pLDDT": True,
    "tcr_1_cdr_2_5_mean_pLDDT": True,
    "tcr_1_cdr_3_mean_pLDDT": True,
    "tcr_2_cdr_1_mean_pLDDT": True,
    "tcr_2_cdr_2_mean_pLDDT": True,
    "tcr_2_cdr_2_5_mean_pLDDT": True,
    "tcr_2_cdr_3_mean_pLDDT": True,
    "mhc_helices_mean_pLDDT": True,
    "tcr_mhc_contacts": True,
    "peptide_tcr_contacts": True,
}
