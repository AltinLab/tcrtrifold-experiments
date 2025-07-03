from .utils import FORMAT_ANTIGEN_COLS


import numpy as np
import polars as pl
from sklearn.metrics import auc
import sklearn.metrics as metrics
from scipy.stats import pearsonr, spearmanr


def train_model(train_df, featnames, model_class, **kwargs):

    n_feat = len(featnames)
    output = ["cognate"]
    train_np = train_df.select(featnames + output).to_numpy()
    X_train = train_np[:, :n_feat]
    y_train = train_np[:, n_feat]
    model = model_class(**kwargs)
    model.fit(X_train, y_train)

    return model


def test_model(model, featnames, test_df):
    n_feat = len(featnames)
    output = ["cognate"]

    test_np = test_df.select(featnames + output).to_numpy()
    X_test = test_np[:, :n_feat]
    y_test = test_np[:, n_feat]

    y_pred_proba = model.predict_proba(X_test)[:, 1]

    fpr, tpr, threshold = metrics.roc_curve(y_test, y_pred_proba)
    roc_auc = metrics.auc(fpr, tpr)
    return fpr, tpr, threshold, roc_auc


def merge_roc_curves(fprs, tprs):
    """
    Merges multiple ROC curves by averaging TPRs at common FPR values.

    Parameters:
    fprs (list of np.array): List of false positive rate arrays.
    tprs (list of np.array): List of true positive rate arrays.

    Returns:
    tuple: Average FPR and TPR arrays.
    """
    # Create a set of all unique FPR values
    all_fpr = np.sort(np.unique(np.concatenate(fprs)))

    # Interpolate TPR values at the unique FPR points
    mean_tprs = np.zeros_like(all_fpr, dtype=np.float64)
    for i, fpr in enumerate(fprs):
        tpr = tprs[i]
        mean_tprs += np.interp(all_fpr, fpr, tpr)

    # Average TPRs and calculate AUC
    mean_tprs /= len(fprs)
    mean_auc = auc(all_fpr, mean_tprs)

    return all_fpr, mean_tprs, mean_auc


def within_antigen_auc(
    triad_df,
    antigen_df,
    featnames,
    model_class,
    model_kwargs,
    n_iter=10,
):
    out_df = []

    for row in antigen_df.iter_rows(named=True):
        antigen = pl.DataFrame([row]).select(pl.exclude("job_name"))
        fprs = []
        tprs = []

        focal_antigen_triads = triad_df.join(antigen, on=FORMAT_ANTIGEN_COLS)

        cognate = focal_antigen_triads.filter(pl.col("cognate")).select(
            featnames + ["cognate"]
        )

        non_cognate = focal_antigen_triads.filter(~pl.col("cognate")).select(
            featnames + ["cognate"]
        )

        # maintain 10:1 neg:pos ratio in sampling
        cognate_test_size = int(cognate.height * 0.1)

        non_cognate_test_size = int(non_cognate.height * 0.1)

        for i in range(n_iter):

            cognate_shuffle = cognate.sample(fraction=1, shuffle=True)
            non_cognate_shuffle = non_cognate.sample(fraction=1, shuffle=True)

            cognate_test, cognate_train = cognate_shuffle.head(
                cognate_test_size
            ), cognate_shuffle.tail(-cognate_test_size)

            non_cognate_test, non_cognate_train = non_cognate_shuffle.head(
                non_cognate_test_size
            ), non_cognate_shuffle.tail(-non_cognate_test_size)

            train = pl.concat([cognate_train, non_cognate_train])
            test = pl.concat([cognate_test, non_cognate_test])

            model = train_model(train, featnames, model_class, **model_kwargs)

            fpr, tpr, threshold, roc_auc = test_model(
                model,
                featnames,
                test,
            )

            fprs.append(fpr)
            tprs.append(tpr)

        fpr, tpr, roc_auc = merge_roc_curves(fprs, tprs)

        row["fpr"] = fpr.tolist()
        row["tpr"] = tpr.tolist()
        row["roc_auc"] = roc_auc

        out_df.append(row)

    return pl.DataFrame(out_df)


def antigen_cross_validation_auc(
    triad_df,
    antigen_df,
    featnames,
    model_class,
    model_kwargs,
    balanced_ratio=False,
):

    out_df = []

    for row in antigen_df.iter_rows(named=True):
        antigen = pl.DataFrame([row]).select(pl.exclude("job_name"))
        focal_antigen_triads = triad_df.join(antigen, on=FORMAT_ANTIGEN_COLS)
        non_focal_antigen_triads = triad_df.join(
            antigen, on=FORMAT_ANTIGEN_COLS, how="anti"
        )

        test, train = focal_antigen_triads, non_focal_antigen_triads

        if balanced_ratio:
            train = pl.concat(
                [
                    train.filter(pl.col("cognate")),
                    train.filter(~pl.col("cognate"))
                    .sample(fraction=1, shuffle=True)
                    .head(train.filter(pl.col("cognate")).height),
                ]
            )

        model = train_model(train, featnames, model_class, **model_kwargs)

        fpr, tpr, threshold, roc_auc = test_model(
            model,
            featnames,
            test,
        )

        new_row = row.copy()

        new_row["fpr"] = fpr.tolist()
        new_row["tpr"] = tpr.tolist()
        new_row["roc_auc"] = roc_auc

        out_df.append(new_row)

    return pl.DataFrame(out_df)
