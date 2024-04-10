#!/usr/bin/python
"""Class for handling correlations in scmodule."""
# --------------------------------------------
# Class for handling correlations in scmodule:
# Updated: 07/02/21
# --------------------------------------------
from .utils_correlation import (
    calculate_correlation,
    calculate_correlation_estimate,
    calculate_correlation_estimate_sd
)

import logging
import numpy as np
import pandas as pd
from scipy import sparse

# TODO: DELETE THIS, WAS GUTTED FOR PARTS

# Standalone correlation handler:
class correlation(object):
    """
    Handles correlation computation on raw + decorrelated data.

    Arguments:
        X: Matrix (n x d) for which to calculate the d x d correlation
        U, s, V: SVD of X matrix (not required).
            U is (n x k), s is (k), and V is in (k x d)
        k: Number of SVD components to compute (default = 100)
        calc_raw: Whether to also return correlation on untransformed data

    Functions:
        get_correlation()
            Obtain the overall correlation matrices

        get_correlation_subset(indices)
            Obtain the correlation matrix for a subset of SVD components
    """

    def __init__(self, X, margin, U, s, V,
                 k=100, calc_raw=False, center=False, mean=None):
        """Initialize correlation handling object."""
        # Input, required (calculated PCA outside of this)
        self.X = X
        self.margin = margin
        self.U = U
        self.s = s
        self.V = V

        # Arguments:
        self.k = k
        self.calc_raw = calc_raw
        self.center = center
        self.mean = mean
        # Output:
        self.corr_est = None
        self.corr_raw = None
        self.use_v = True
        logging.debug("X: " + str(self.X.shape))
        logging.debug("Margin: " + str(self.margin.shape))


    def get_correlation(self, indices=None, keep_first=False, power=0):
        """Get correlation estimates from SVD and raw if calc_raw is True."""

        if indices is None:


        else:
            corr_subset = calculate_correlation_estimate(
                self.U, self.s, self.V, power=power, indices=indices)


        if self.corr_est is None or \
                self.power != power:
            self._calculate_estimated_correlation(
                keep_first=keep_first, power=power)

        self.corr_raw = self.get_raw_correlation()

        return (self.corr_est, self.corr_raw)

    def get_raw_correlation(self, force=False):
        if force:
            self.calc_raw = True
        if self.corr_raw is None and self.calc_raw:
            self.corr_raw = calculate_correlation(self.X, center=self.center)
        # NOTE: else corr_raw is already None
        return(self.corr_raw)

#     def get_correlation_subset(self, indices, power=0):
#         """Get a correlation estimate from a subset of SVD columns."""
#         logging.info("Calculating correlation from subset of SVD columns")
#         corr_subset = calculate_correlation_estimate(
#                 self.U, self.s, self.V, power=power, indices=indices)
#         return corr_subset

    def _calculate_estimated_correlation(self, power=0, keep_first=False):
        self.power = power
        if self.use_v:
            # Keep or remove first component (defaults to removing):
            indices = None if keep_first else np.arange(1, self.U.shape[1])
            logging.debug(indices)
            self.corr_est = calculate_correlation_estimate(
                self.U, self.s, self.V, power=power, indices=indices)

    # TODO: allow both types of adjustments (centralize in graph vs. not
    # in graph types)

    # TODO: Estimate the raw correlation SD as well.
    def estimate_corr_sd(self, nperm=25, seed=1):
        corr_mean, corr_sd = calculate_correlation_estimate_sd(
            self.U, self.s, self.V, nperm=nperm, seed=seed)
        return (corr_mean, corr_sd)
