import matplotlib.pyplot as plt
from torch_geometric.data import Data


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
