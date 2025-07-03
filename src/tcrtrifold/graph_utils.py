import matplotlib.pyplot as plt


def plot_auc_per_antigen(
    cv_df, ax=None, title=None, id_cols=[], roc_name="roc_auc"
):
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
