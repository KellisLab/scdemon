#!/usr/bin/env python3
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


def calculate_svd_covar_corr(ut, obsdf, cvlist,
                             covariate_matrices={}, force_run=False):
    """Calculate covariate correlations with SVD."""
    for covar in cvlist:
        if covar not in covariate_matrices.keys() or force_run:
            covariate_col = obsdf[covar]
            covariate_matrices[covar] = calculate_svd_covar_corr_single(
                ut, covariate_col)
    return(covariate_matrices)


def calculate_svd_covar_corr_single(ut, covariate_col):
    if covariate_col.dtype.name == "category":
        covar_dummy = pd.get_dummies(covariate_col)
        corr = np.zeros((ut.shape[0], covar_dummy.shape[1]))
        for i, covariate_lvl in enumerate(covar_dummy.columns):
            corr[:, i] = vcorrcoef(
                ut, covar_dummy[covariate_lvl].to_numpy().T)
    else:
        covariate_col = covariate_col.to_numpy()
        if np.sum(np.isnan(covariate_col)) > 0:
            ind = np.where((1 - np.isnan(covariate_col)) == 1)[0]
            corr = vcorrcoef(
                ut[:, ind], covariate_col[ind, np.newaxis].T)[:, np.newaxis]
        else:
            corr = vcorrcoef(ut, covariate_col[:, np.newaxis].T)[:, np.newaxis]
    return(corr)

