#!/usr/bin/python
"""Class for handling an adjacency matrix."""
from .utils_adjacency import (
    zscore_from_bivariate_cutoff,
    adj_from_sd_est
)
from .utils_pruning import (
    prune_adj_degree,
    prune_adj_scale,
    prune_adj_knn,
)


import logging
import numpy as np
from scipy import sparse


class adjacency_matrix(object):
    """
    Turn a correlation matrix into an adjacency matrix.

    Process given correlation matrix into the adjacency matrix.
    """

    def __init__(
        self,
        corr,
        adjacency=None,
        method='bivariate',
        corr_sd=None,
        labels=None,
        margin=None,
        # Threshold arguments:
        cutoff=0.4, z=4.5,
        zero_outliers=True,
        keep_all_z=False,
        # Pruning graph arguments:
        knn_k=None,  # k-NN
        scale=None,
        degree_cutoff=0
    ):
        """Initialize adjacency matrix class."""
        self.corr = corr
        self.adjacency = adjacency
        self.method = method
        self.corr_sd = corr_sd  # Estimated sd for each pair
        self.labels = labels
        self.indices = np.arange(len(self.labels))
        self.margin = margin  # Margin (fraction non-zero) of dataset

        # Arguments for thresholding
        self.cutoff = cutoff
        self.z = z  # For adjacency cutoff
        self.zero_outliers = zero_outliers  # For getting estimates of cutoffs
        self.keep_all_z = keep_all_z

        # Arguments for pruning
        self.knn_k = knn_k
        self.scale = scale
        self.degree_cutoff = degree_cutoff
        logging.debug("Margin: " + str(self.margin.shape))
        logging.debug("Corr: " + str(self.corr.shape))

        # Create adjacency matrix directly with constructor:
        # Also prunes adjacency, run once.
        self._create_adjacency_matrix()

    def _create_adjacency_matrix(self):
        """Threshold + prune correlation to create adjacency matrix."""
        # Threshold the correlation matrix to get the adjacency
        if self.adjacency is None:
            self.adjacency = self._threshold_correlation_matrix()

        # Full and mostly empty adjacency, for multiplexing graphs:
        self.full_adjacency = self.adjacency
        self.full_labels = self.labels

        # Prune adjacency matrix:
        self._prune_adjacency()

    def get_adjacency(self):
        """Get the adjacency matrix and the final kept labels."""
        return(self.adjacency, self.labels, self.indices)

    def _threshold_correlation_matrix(self):
        """Threshold correlation to get adjacency matrix."""
        logging.info(
            "Thresholding correlation with method '%s' to make adjacency" %
            self.method)
        if self.method == 'cutoff':
            adjacency = self.corr * (self.corr > self.cutoff)
            adjacency = adjacency - np.diag(np.diag(adjacency))
            adjacency = sparse.coo_matrix(adjacency)
        elif self.method == 'bivariate':
            _, adjacency, _ = \
                zscore_from_bivariate_cutoff(
                    self.corr, self.margin, z=self.z,
                    zero_outliers=self.zero_outliers,
                    keep_all_z=self.keep_all_z)
        elif self.method == 'sd':
            adjacency = adj_from_sd_est(
                self.corr, self.corr_sd, self.margin, z=self.z)
        else:
            raise ValueError(
                "Method %s not in (cutoff, bivariate, sd)" %
                self.method)
        return(adjacency)

    def _prune_adjacency(self):
        # Prune the created matrix:
        if self.scale is not None:
            self._prune_by_scale()
        if self.knn_k is not None:
            self._prune_by_knn()
        # Clean up nodes with low or 0 degree after pruning by other methods:
        self._prune_by_degree()

    def _prune_by_scale(self):
        """Prune adjacency by a relative cutoff of outgoing edges."""
        logging.info("Pruning by relative cutoff, scale=" + str(self.scale))
        self.adjacency = prune_adj_scale(self.adjacency, scale=self.scale)

    def _prune_by_knn(self):
        """Prune adjacency by a k-NN keeping top k edges."""
        logging.info("Pruning by KNN, K=" + str(self.knn_k))
        self.adjacency = prune_adj_knn(self.adjacency, k=self.knn_k)

    def _prune_by_degree(self):
        """Prune adjacency by the degree of each node."""
        logging.info("Pruning by degree, cutoff=" + str(self.degree_cutoff))
        # Remove nodes with no/very few links:
        self.adjacency, ind = prune_adj_degree(
            self.adjacency, self.degree_cutoff)
        self.labels = self.labels[ind]
        self.indices = self.indices[ind]
        self.adjacency = sparse.csr_matrix(self.adjacency)
