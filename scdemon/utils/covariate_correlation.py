#!/usr/bin/python
"""Utility functions for correlation between covariates and PCs"""
import logging
import numpy as np
import pandas as pd


def vcorrcoef(X, y):
    """Vectorized correlation, matrix vs. vector."""
    Xm = np.reshape(np.mean(X, axis=1), (X.shape[0], 1))
    ym = np.mean(y)
    r_num = np.sum(np.array(X - Xm) * (y - ym), axis=1)
    r_den = np.sqrt(np.sum(np.array(X - Xm) ** 2, axis=1)
                    * np.sum((y - ym) ** 2))
    r = r_num / r_den
    return r


def calculate_svd_covar_corr(ut, obsdf, cvlist, cv_mats={}, force_run=False):
    """Calculate covariate correlations with SVD."""
    for covar in cvlist:
        if covar not in cv_mats.keys() or force_run:
            cvcol = obsdf[covar]
            cv_mats[covar] = calculate_svd_covar_corr_single(ut, cvcol)
    return(cv_mats)


def calculate_svd_covar_corr_single(ut, cvcol):
    if cvcol.dtype.name == "category":
        covar_dummy = pd.get_dummies(cvcol)
        corr = np.zeros((ut.shape[0], covar_dummy.shape[1]))
        for i, cv_lvl in enumerate(covar_dummy.columns):
            corr[:, i] = vcorrcoef(ut, covar_dummy[cv_lvl].to_numpy().T)
    else:
        cvcol = cvcol.to_numpy()
        if np.sum(np.isnan(cvcol)) > 0:
            ind = np.where((1 - np.isnan(cvcol)) == 1)[0]
            corr = vcorrcoef(
                ut[:, ind], cvcol[ind, np.newaxis].T
            )[:, np.newaxis]
        else:
            corr = vcorrcoef(ut, cvcol[:, np.newaxis].T)[:, np.newaxis]
    return(corr)

