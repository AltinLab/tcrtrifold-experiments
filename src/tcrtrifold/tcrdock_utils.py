import polars as pl
import numpy as np
from scipy.stats import chi2


def cossin_embed(X):
    X_new = np.empty((X.shape[0], 7))
    X_new[:, 0] = X[:, 0]
    X_new[:, 1] = np.sin(X[:, 1])
    X_new[:, 2] = np.cos(X[:, 1])
    X_new[:, 3:] = X[:, 2:]

    return X_new


def un_cossin_embed(X):
    X_new = np.empty((X.shape[0], 6))
    X_new[:, 0] = X[:, 0]
    X_new[:, 1] = np.arctan2(X[:, 1], X[:, 2])
    X_new[:, 2:] = X[:, 3:]

    return X_new


def dgeom_ndarr_from_dgeom_series(dgeom_series):
    dgeom_ndarr = np.array(
        pl.DataFrame({"dgeom": dgeom_series})
        .with_columns(
            pl.concat_list(
                [
                    pl.col("dgeom").struct.field("d"),
                    pl.col("dgeom").struct.field("torsion"),
                    pl.col("dgeom").struct.field("tcr_unit_y"),
                    pl.col("dgeom").struct.field("tcr_unit_z"),
                    pl.col("dgeom").struct.field("mhc_unit_y"),
                    pl.col("dgeom").struct.field("mhc_unit_z"),
                ]
            ).alias("dgeom")
        )
        .select(pl.col("dgeom").implode())
        .to_series()
        .to_list()[0]
    )
    return cossin_embed(dgeom_ndarr)


def mn_distr_from_dgeom_ndarr(dgeom_ndarr):
    # dgeom_ndarr = dgeom_ndarr_from_dgeom_series(dgeom_ndarr)

    mu = dgeom_ndarr.mean(axis=0)
    cov = np.cov(dgeom_ndarr, rowvar=False)
    # cov += np.eye(cov.shape[0]) * 1e-6
    invcov = np.linalg.inv(cov)

    return mu, invcov


def mn_distance_from(dgeom_ndarr, mu, invcov):
    # dgeom_ndarr = dgeom_ndarr_from_dgeom_series(dgeom_series)

    diff = dgeom_ndarr - mu
    dm2 = np.sum((diff @ invcov) * diff, axis=1)
    p = chi2.sf(dm2, df=7)

    return np.sqrt(dm2), p
