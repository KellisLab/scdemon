#!/usr/bin/env python3
"""Utility functions for calculating / estimating correlation."""
import logging
import numpy as np


def triu_mask(A):
    """Get the upper triangular part a matrix."""
    m = A.shape[0]
    r = np.arange(m)
    mask = r[:, None] < r
    return A[mask]


def calculate_correlation(X, center=False):
    """Calculate the correlation of columns of a matrix."""
    logging.info("Calculate correlation of columns of matrix X")
    # TODO: CENTER (+ OTHER WAYS of calc)
    # TODO: Allow numpy corrcoeff
    if center:
        # TODO: See difference with centering X:
        raise NotImplementedError("Raw centering not implemented")
        corr = None
    else:
        from scipy.sparse import issparse
        xTx = X.T.dot(X) / X.shape[0]
        if issparse(xTx):
            xTx = xTx.toarray()
        sd = np.sqrt(np.diag(xTx))
        xTx = xTx / sd[:, np.newaxis]
        corr = xTx / sd[np.newaxis, :]
        del(xTx, sd)
    return(corr)


def calculate_correlation_estimate(U, s, V, power=0, indices=None):
    if type(power) is list:
        corr_est = np.zeros((V.shape[1], V.shape[1]))
        j = 0
        for pw in power:
            j += 1
            corr_est += calculate_single_correlation_estimate(
                U=U, s=s, V=V, power=pw, indices=indices)
        corr_est = corr_est / (j * 1.0)
    else:
        corr_est = calculate_single_correlation_estimate(
            U=U, s=s, V=V, power=power, indices=indices)
    return(corr_est)


def calculate_single_correlation_estimate(U, s, V, power=0, indices=None):
    """Calculate an SVD-derived estimate of the correlation matrix."""
    logging.info("Estimating correlation of columns of matrix X with its SVD")
    scale = U.shape[0] * V.shape[1]
    if indices is not None:
        logging.debug("Reducing to subset of indices")
        s = s[indices]
        V = V[indices, :]
    if power != 0:
        logging.debug(f"Using power {power}")
        s_red = np.diag(s**power)
        X_cov = V.T.dot(s_red).dot(V) / scale
    else:
        X_cov = V.T.dot(V) / scale

    X_sd = np.sqrt(np.diag(X_cov))
    cv = X_cov / X_sd[:, np.newaxis]
    corr_est = cv / X_sd[np.newaxis, :]
    del (cv, X_cov, X_sd, scale)
    return(corr_est)


def calculate_correlation_estimate_sd(U, s, V, nperm=25, seed=1):
    """Estimate the sd of correlations by subsampling V."""
    np.random.seed(seed)  # Reproducible sampling
    k = U.shape[1]
    full_ind = np.random.randint(k, size=(k // 2, nperm))
    corr_est = np.zeros((V.shape[1], V.shape[1], nperm))
    for i in range(nperm):
        logging.debug(i)
        corr_est[:, :, i] = calculate_correlation_estimate(
            U, s, V, indices=full_ind[:, i])

    corr_mean = np.mean(corr_est, axis=-1)
    corr_sd = np.std(corr_est, axis=-1)
    del(corr_est, full_ind)
    return (corr_mean, corr_sd)

