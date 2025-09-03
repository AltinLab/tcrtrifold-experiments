import matplotlib.pyplot as plt
from torch_geometric.data import Data
import numpy as np
import polars as pl
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import ttest_rel, ttest_1samp


def plot_tg_data_3d(graph: Data):
    node_color = "red"
    edge_color = "blue"

    pos = graph.pos.cpu().numpy()
    edge_index = graph.edge_index.cpu().numpy()
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")
    xs, ys, zs = pos[:, 0], pos[:, 1], pos[:, 2]

    color_int = [ord(segid) - 65 for segid in graph.segids]
    ax.scatter(xs, ys, zs, s=50, c=color_int, cmap="prism")

    for src, dst in edge_index.T:
        x_line = [pos[src, 0], pos[dst, 0]]
        y_line = [pos[src, 1], pos[dst, 1]]
        z_line = [pos[src, 2], pos[dst, 2]]
        ax.plot(x_line, y_line, z_line, c=edge_color, alpha=0.5, linewidth=0.5)

    plt.tight_layout()
    plt.show()


def plot_auc_per_antigen(cv_df, ax=None, title=None, id_cols=[], roc_name="roc_auc"):
    """
    Plot ROC curves from antigen_cross_validation_auc output.

    Parameters
    ----------
    cv_df : pl.DataFrame
        Must have List-type columns 'fpr' and 'tpr', and a float column 'roc_auc',
        plus one or more identifier columns (e.g. peptide, mhc_1, mhc_2).
    ax : matplotlib.axes.Axes, optional
        If provided, plot into this Axes. Otherwise creates a new figure+axis.
    title : str, optional
        Overall title for the plot.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The Axes with your ROC curves.
    """
    # create new figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    aucs = []

    cv_df = cv_df.sort(by=id_cols)

    # plot each ROC
    for row in cv_df.iter_rows(named=True):
        fpr = row["fpr"]
        tpr = row["tpr"]
        auc = row[roc_name]
        aucs.append(auc)
        # build label from the ID columns
        label_parts = [f"{col}: {row[col]}" for col in id_cols]
        label = ", ".join(label_parts) + f" (AUC {auc:.2f})"
        ax.plot(fpr, tpr, lw=1.5, label=label)

    mean_auc = sum(aucs) / len(aucs)
    ax.text(
        0.05,
        0.95,
        f"Mean AUC: {mean_auc:.2f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize="small",
        bbox=dict(boxstyle="round,pad=0.3", alpha=0.3),
    )

    # chance line
    ax.plot([0, 1], [0, 1], "--", lw=1, color="grey")

    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    if title:
        ax.set_title(title)
    ax.legend(loc="right", fontsize=5)
    ax.grid(True)

    return ax


def per_feat_auc(df, featnames_dict, title, axes=None, col=None):
    """
    Plot per-feature AUC boxplots, optionally into provided axes grid.

    Parameters
    ----------
    df : pl.DataFrame
        DataFrame containing a boolean 'cognate' column and feature columns.
    featnames_dict : dict
        Mapping of feature names to a bool flag (unused here but preserved for API).
    title : str
        Title for the set of plots.
    axes : array-like of Axes or 2D array of Axes, optional
        If provided, draws plots into these axes instead of creating new ones.
    col : int, optional
        If `axes` is a 2D array, select this column (0-indexed) for plotting.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing the plots.
    """
    featnames = list(featnames_dict.keys())

    # Determine axes and figure
    if axes is None:
        fig, ax_arr = plt.subplots(nrows=len(featnames), ncols=1, figsize=(20, 80))
    else:
        # axes provided by user
        fig = axes[0].figure if not hasattr(axes, "shape") else axes[0, 0].figure
        if hasattr(axes, "shape") and axes.ndim == 2:
            if col is None:
                raise ValueError("Must specify 'col' when passing a 2D axes array.")
            ax_arr = axes[:, col]
        else:
            ax_arr = axes

    # Plot each feature
    for ax, feat in zip(ax_arr, featnames):
        noncognate = df.filter(~pl.col("cognate")).select(feat).to_series().to_numpy()
        cognate = df.filter(pl.col("cognate")).select(feat).to_series().to_numpy()

        # stats
        t_stat, p_value = ttest_ind(
            cognate, noncognate, equal_var=False, alternative="two-sided"
        )
        auc = roc_auc_score(
            [1] * len(cognate) + [0] * len(noncognate),
            list(cognate) + list(noncognate),
        )

        # boxplot
        bp = ax.boxplot(
            [noncognate, cognate],
            positions=[2, 1],
            vert=False,
            patch_artist=True,
            widths=0.4,
        )
        colors = ["red", "green"]
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
        for median in bp["medians"]:
            median.set(color="black", linewidth=3)

        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xlabel(f"p = {p_value:.2g}\nAUC = {auc:.2f}", size="large")

    # Title and labels
    ax_arr[0].set_title(title, fontsize="large")

    if col is None or col == 0:
        for ax, row in zip(ax_arr, featnames):
            ax.set_ylabel(row, rotation=0, size="large", labelpad=100)

        handles = [
            mpatches.Patch(facecolor="red", label="Noncognate"),
            mpatches.Patch(facecolor="green", label="Cognate"),
        ]
        fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.15, 0.95))

    plt.tight_layout()
    return fig


def heatmap(df, suptitle=None):

    contact_maps_cognate = np.array(
        df.filter(pl.col("cognate")).select("contact_map").to_series().to_list()
    )

    contact_maps_noncognate = np.array(
        df.filter(~pl.col("cognate")).select("contact_map").to_series().to_list()
    )
    cmap_gray_blue = LinearSegmentedColormap.from_list(
        "WhiteToBlue", ["white", "darkblue"]
    )

    mean_cognate = np.mean(contact_maps_cognate, axis=0)
    mean_noncognate = np.mean(contact_maps_noncognate, axis=0)

    col_labels = [str(i) for i in range(1, 10)] + ["HLA A", "HLA B"]

    row_segments = []
    for chain in ["Alpha", "Beta"]:
        for seg in [
            "fwr_1",
            "cdr_1",
            "fwr_2",
            "cdr_2",
            "fwr_3",
            "cdr_3",
            "fwr_4",
        ]:
            row_segments.append(f"{' '.join(seg.split('_')).upper()}")

    fig, axs = plt.subplots(2, 1, figsize=(8, 10), constrained_layout=True)

    vmin = min(mean_cognate.min(), mean_noncognate.min())
    vmax = max(mean_cognate.max(), mean_noncognate.max())

    t = []
    p = []
    log = []

    if df.select("mhc_class")[0].item() == "II":
        exp_pos = [1, 2, 4, 6, 7]
        non_exp_pos = [0, 3, 5, 8]

    else:
        exp_pos = [3, 4, 5, 6, 7]
        non_exp_pos = [0, 1, 2, 8]

    for ax, data, title, full_data in zip(
        axs,
        [mean_cognate, mean_noncognate],
        ["Cognate", "Noncognate"],
        [contact_maps_cognate, contact_maps_noncognate],
    ):
        im = ax.imshow(data, aspect="auto", cmap=cmap_gray_blue, vmin=vmin, vmax=vmax)
        ax.set_title(title + f" (n={full_data.shape[0]})", pad=20)
        ax.set_xticks(np.arange(len(col_labels)))
        ax.set_xticklabels(col_labels)
        ax.set_yticks(np.arange(len(row_segments)))
        ax.set_yticklabels(row_segments)

        # ax.set_xticks(range(mean_cognate.shape[1]), labels=col_labels,
        #               ha="right", rotation_mode="anchor")
        # ax.set_yticks(range(mean_cognate.shape[0]), labels=row_segments)

        ax.spines[:].set_visible(False)

        ax.grid(which="minor", color="w", linestyle="-", linewidth=3)
        ax.tick_params(which="minor", bottom=False, left=False)

        # ax.tick_params(axis='both', length=0)

        ax2 = ax.twinx()
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks([3, 10])
        ax2.set_yticklabels(["alpha", "beta"])
        ax2.spines["left"].set_position(("outward", 60))
        ax2.spines["left"].set_visible(False)
        ax2.spines["right"].set_visible(False)
        ax2.yaxis.set_ticks_position("left")
        ax2.yaxis.set_label_position("left")
        ax2.tick_params(axis="y", length=0)

        # full data = n_samples, 14, 11

        exp_arr = full_data[:, :, exp_pos].sum(axis=1).sum(axis=1) / len(exp_pos)
        non_exp_arr = full_data[:, :, non_exp_pos].sum(axis=1).sum(axis=1) / len(
            non_exp_pos
        )

        t_stat, p_value = ttest_rel(
            exp_arr,
            non_exp_arr,
            # equal_var=False,
            alternative="greater",
        )
        # # tcr = np.sum(full_data[:, [5, 12]][:, :, :9], axis=1)
        # # tcr_mean = np.mean(tcr, axis=0)
        logfold = np.log2(exp_arr.mean() / non_exp_arr.mean())

        # log_exp = np.log2(exp_arr + 1e-6)
        # log_non = np.log2(non_exp_arr + 1e-6)

        # t_stat, p_value = ttest_rel(log_exp, log_non, alternative="greater")
        # log2_fold_change = (log_exp - log_non).mean()

        t.append(t_stat)
        p.append(p_value)
        log.append(logfold)

    print(
        f"Cognate Log2FC in TCR-facing positions ({[i + 1 for i in exp_pos]}) vs MHC-facing posititions ({[i + 1 for i in non_exp_pos]}): {log[0]:.2f}\np-value: {p[0]:.2e}",
    )
    print(
        f"Non-cognate Log2FC in TCR-facing positions: {log[1]:.2f}\np-value: {p[1]:.2e}",
    )

    # fig.text(
    #     1.03,
    #     0.9,
    #     f"Log2FC in TCR-facing positions: {log[0]:.2f}\np-value: {p[0]:.2e}",
    #     fontsize="large",
    #     bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="black"),
    # )
    # fig.text(
    #     1.03,
    #     0.4,
    #     f"Log2FC in TCR-facing positions: {log[1]:.2f}\np-value: {p[1]:.2e}",
    #     fontsize="large",
    #     bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="black"),
    # )

    cbar = fig.colorbar(im, ax=axs, orientation="vertical", pad=0.02)
    cbar.set_label("Mean num contacts", rotation=-90, va="bottom")

    if suptitle:
        fig.suptitle(suptitle, fontsize="large", y=1.03)

    return fig


def rotation_matrix_from_vectors(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Compute the rotation matrix that aligns vector a to vector b
    using Rodrigues' rotation formula (proper rotation, no reflection).
    """
    a_unit = a / np.linalg.norm(a)
    b_unit = b / np.linalg.norm(b)
    v = np.cross(a_unit, b_unit)
    s = np.linalg.norm(v)
    c = np.dot(a_unit, b_unit)
    # skew-symmetric cross-product matrix
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    if s < 1e-8:
        # vectors are parallel or antiparallel
        if c > 0:
            return np.eye(3)
        else:
            # 180° rotation: choose any orthogonal axis
            # e.g. perpendicular to a_unit
            ortho = np.array([1, 0, 0])
            if abs(a_unit @ ortho) > 0.9:
                ortho = np.array([0, 1, 0])
            axis = np.cross(a_unit, ortho)
            axis /= np.linalg.norm(axis)
            K2 = np.array(
                [
                    [0, -axis[2], axis[1]],
                    [axis[2], 0, -axis[0]],
                    [-axis[1], axis[0], 0],
                ]
            )
            return np.eye(3) + 2 * K2 @ K2
    # Rodrigues formula
    return np.eye(3) + K + (K @ K) * ((1 - c) / (s**2))


def rodrigues_axis_angle(axis: np.ndarray, theta: float) -> np.ndarray:
    """
    Rotation matrix for a rotation of angle theta about a given axis.
    """
    k = axis / np.linalg.norm(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    c, s = np.cos(theta), np.sin(theta)
    return np.eye(3) + s * K + (1 - c) * (K @ K)


def signed_dihedral(u_axis, v1, v2):
    """
    Compute the signed angle from v1→v2 around the unit-axis u_axis
    using atan2(cross, dot).
    All inputs are 3‐vectors.
    """

    # project v1, v2 onto plane ⟂ u_axis
    def proj_plane(v, n):
        return v - n * (n.dot(v))

    p1 = proj_plane(v1, u_axis)
    p2 = proj_plane(v2, u_axis)
    p1 /= np.linalg.norm(p1)
    p2 /= np.linalg.norm(p2)
    # signed angle
    x = np.dot(p1, p2)
    y = np.dot(u_axis, np.cross(p1, p2))
    return np.arctan2(y, x)


def build_tcr_orientation_rotation(
    torsion: float,
    tcr_unit_y: float,
    tcr_unit_z: float,
    mhc_unit_y: float,
    mhc_unit_z: float,
) -> np.ndarray:
    """
    Given 6 scalars:
      - torsion angle (radians) between MHC-y and TCR-z axes
      - (tcr_unit_y, tcr_unit_z): the y,z components of the TCR→MHC direction
         in MHC coordinates
      - (mhc_unit_y, mhc_unit_z): the y,z components of the MHC→TCR direction
         in TCR coordinates

    Returns the 3×3 rotation matrix R_total that:
      1) Aligns the MHC→TCR vector in the TCR frame to the TCR→MHC vector
         in the MHC frame (first docking alignment).
      2) Applies the specified torsion (twist) around that docking axis.

    After this, R_total can be applied to the identity basis:
        rotated_axes = R_total @ np.eye(3)
    whose columns are the TCR x,y,z axes expressed in the MHC frame, with
    the correct docking direction and torsion.
    """
    # 1) Recover full 3D unit docking vectors
    #    TCR in MHC frame:
    ty, tz = tcr_unit_y, tcr_unit_z
    tx = np.sqrt(max(0.0, 1 - ty**2 - tz**2))
    v = np.array([tx, ty, tz])  # unit vector MHC→TCR in MHC coords

    #    MHC in TCR frame:
    my, mz = mhc_unit_y, mhc_unit_z
    mx = np.sqrt(max(0.0, 1 - my**2 - mz**2))
    u = np.array([mx, my, mz])  # unit vector TCR→MHC in TCR coords

    R1 = rotation_matrix_from_vectors(u, -v)
    axis = v

    mhc_y = np.array([0, 1, 0])
    tcr_z = np.array([0, 0, 1])
    z1 = R1 @ tcr_z
    current = signed_dihedral(
        axis,
        mhc_y,
        z1,
    )

    twist = torsion - current

    # 6) build that twist‐rotation about the docking axis
    R2 = rodrigues_axis_angle(axis, twist)
    # 4) Combined rotation
    R_total = R2 @ R1
    # R_total = R1
    return (R_total, v, u, axis)


def plot_docking_geometry(
    d, torsion, tcr_unit_y, tcr_unit_z, mhc_unit_y, mhc_unit_z, ax=None
):
    """
    3D visualization of docking geometry:
      - MHC at origin, with its axes drawn (x up, y right, z forward)
      - TCR at distance d in the MHC frame, direction by (tcr_unit_y, tcr_unit_z)
      - A dashed line connects MHC→TCR
      - An arc in the MHC yz-plane shows the torsion angle
      - TCR with its axes draw

      TCR's axes are determined by aligning the vector describing the position of the
      TCR from the MHC's frame and the vector describing the position of the MHC from
      the TCR's frame

    """
    if ax is None:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection="3d")
        first_plot = True
    else:
        fig = ax.get_figure()
        first_plot = False

    tcr_axes, tcr_unit, mhc_unit, bond_axis = build_tcr_orientation_rotation(
        torsion, tcr_unit_y, tcr_unit_z, mhc_unit_y, mhc_unit_z
    )

    tcr_pos = np.array(tcr_unit) * d
    tcr_axes = tcr_axes.T

    axis_length = d * 0.25

    if first_plot:
        ax.scatter(0, 0, 0, color="blue", s=50, label="MHC (origin)")
        mhc_axes = [
            ("red", np.array([1, 0, 0])),
            ("green", np.array([0, 1, 0])),
            ("blue", np.array([0, 0, 1])),
        ]
        for color, vec in mhc_axes:
            ax.quiver(0, 0, 0, *(vec * 5), color=color)

    if first_plot:
        label = "TCR"
    else:
        label = None

    ax.scatter(*tcr_pos, color="red", s=50, label=label)

    tcr_axes = [
        ("red", tcr_axes[0]),
        ("green", tcr_axes[1]),
        ("blue", tcr_axes[2]),
    ]

    for color, vec in tcr_axes:
        ax.quiver(*tcr_pos, *(vec * 5), color=color)

    if first_plot:
        label = "bond axis"
    else:
        label = None
    ax.plot(
        [0, tcr_pos[0]],
        [0, tcr_pos[1]],
        [0, tcr_pos[2]],
        "k--",
        label=label,
    )

    # Draw torsion arc in MHC yz-plane
    # angles = np.linspace(0, torsion, 100)
    # y_unitvec = np.array([0, 1, 0])
    # z_unitvec = np.array([0, 0, 1])
    # arc = np.outer(np.cos(angles), y_unitvec) + np.outer(
    #     np.sin(angles), z_unitvec
    # )
    # arc *= axis_length * 0.3
    # ax.plot(arc[:, 0], arc[:, 1], arc[:, 2], color="green", label="torsion")
    if not first_plot:
        curr_ylim = list(ax.get_ylim())
        curr_xlim = list(ax.get_xlim())
        curr_zlim = list(ax.get_zlim())

    # Equalize axis ranges
    all_pts = np.vstack([[0, 0, 0], tcr_pos])
    max_range = (all_pts.max(axis=0) - all_pts.min(axis=0)).max() * 1.1
    mid = (all_pts.max(axis=0) + all_pts.min(axis=0)) / 2

    if first_plot:
        ax.set_xlim(mid[0] - max_range / 2, mid[0] + max_range / 2)
        ax.set_ylim(mid[1] - max_range / 2, mid[1] + max_range / 2)
        ax.set_zlim(mid[2] - max_range / 2, mid[2] + max_range / 2)

    else:
        if curr_ylim[0] > mid[1] - max_range / 2:
            curr_ylim[0] = mid[1] - max_range / 2
        if curr_ylim[1] < mid[1] + max_range / 2:
            curr_ylim[1] = mid[1] + max_range / 2
        if curr_xlim[0] > mid[0] - max_range / 2:
            curr_xlim[0] = mid[0] - max_range / 2
        if curr_xlim[1] < mid[0] + max_range / 2:
            curr_xlim[1] = mid[0] + max_range / 2
        if curr_zlim[0] > mid[2] - max_range / 2:
            curr_zlim[0] = mid[2] - max_range / 2
        if curr_zlim[1] < mid[2] + max_range / 2:
            curr_zlim[1] = mid[2] + max_range / 2

        ax.set_xlim(curr_xlim)
        ax.set_ylim(curr_ylim)
        ax.set_zlim(curr_zlim)

    ax.set_xlabel("MHC x")
    ax.set_ylabel("MHC y")
    ax.set_zlabel("MHC z")
    ax.legend()
    fig.tight_layout()

    # ax.view_init(elev=20, azim=45)
    # ax.view_init(elev=0, azim=180)
    ax.view_init(azim=-45, vertical_axis="x")

    return ax
