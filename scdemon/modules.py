#!/usr/bin/env python3
"""Class for computing modules on adata object."""

# Internal imports:
from .graph import gene_graph, adjacency_matrix
from .utils import (
    # Correlation:
    calculate_correlation,
    calculate_correlation_estimate,
    calculate_correlation_estimate_sd,
    calculate_svd_covar_corr, # Covariates
    # Multigraph
    make_graphlist, partition_graphlist,
)

import os
import re
import gc
import logging

import numpy as np
import pandas as pd
from scipy import sparse
from pandas.api.types import is_categorical_dtype

from anndata import AnnData

class modules(object):
    """
    Calculate gene modules from anndata or stand-alone matrix.
    """

    def __init__(self,
                 adata,
                 X=None, meta=None, # TODO: Add alternative if no adata.obs
                 U=None, s=None, V=None, # If we want to pass these in
                 suffix=None,
                 # Arguments for setup:
                 seed=1,
                 filter_expr=0.05,
                 svd_k=100, # TODO: change to k, change example code to match
                 keep_first_PC=False,
                 process_covariates=False, # Default of processing covars
                 cv_cutoff=0.4,
                 ):
        """
        Initialize gene modules object.
        ---------------

        test
        :param adata: AnnData object containing single-cell dataset

        """
        # TODO: implement stand-alone matrix
        # Dataset:
        self.adata = adata
        self.X = X

        # Reduction
        self.U = U
        self.s = s
        self.V = V

        # Metadata
        self.meta = meta
        self.cvlist = self.adata.obs.columns # TODO: use meta if using X
        self.genes = self.adata.var_names  # labels, TODO different if X

        # Tagline
        self.suffix = suffix if suffix is not None else ''

        # Arguments for setup steps (up to correlation estimate)
        self.k = svd_k
        self.filter_expr = filter_expr
        self.seed = seed
        self.keep_first_PC = keep_first_PC
        self.process_covariates = process_covariates
        self.cv_cutoff = cv_cutoff

        # Additional metadata / files:
        self.graphs = {}  # Store computed graphs
        self.cv_mats = {} # For correlation of covariates vs. PCs
        self.cv_pc_ind = {} # For filtering PCs by covariate

        # Initialize with the random seed
        np.random.seed(self.seed)

    # TODO: Split fixed vs. repeated tasks (for multiple graphs)

    # Part I. Setup data > PCs > correlation estimate
    # -----------------------------------------------
    # TODO: Allow options:
    # TODO: Possibly separate into setup + calculate correlation,
    # so can we can run in multiple times, make multiple graphs
    def setup(self):
        """Set up the dataset and get the correlation estimate."""
        # 1a. normalize (done outside of this)
        self._filter_by_genes() # 1b. subset to genes above expr cutoff
        self._calculate_PCA() # 2. perform SVD or get pre-computed SVD
        if self.process_covariates:
            self._setup_covariate_PC_comparison() # 2a. to filter components
        self._adjust_PCs() # 3. PC adjustment # TODO: Put before or after 2a?

    # NOTE: Filter in place on modules object
    def _filter_by_genes(self):
        """Filter dataset to genes expressed above a given % cutoff."""
        self.margin = calculate_margin_genes(self.adata.X)
        if self.filter_expr is not None:
            above_cutoff = self.margin >= self.filter_expr
            filtgenes = self.adata.var_names[above_cutoff]
            self.margin = self.margin[above_cutoff]
            self.adata = self.adata[:, filtgenes]
            self.genes = self.adata.var_names
            logging.debug("Subset to %s genes." % len(self.genes))

    def _check_PCA(self):
        """Check conditions for re-calculating PCA."""
        if self.U is None or self.s is None or \
                self.V is None or self.U.shape[1] < self.k:
            return False
        else:
            return True

    def _calculate_PCA(self, center=False):
        """
        Calculate PCA components U, s, V if they are not already computed.
        ---
        Defaults to X_pca in adata.obsm if available.
        Otherwise uses sc.tl.pca
        """
        if not self._check_PCA():
            if ("X_pca" not in self.adata.obsm.keys()) or \
                    (self.adata.obsm['X_pca'].shape[1] < self.k):
                # Compute PCA with scanpy options:
                import scanpy as sc
                logging.debug("Computing PCA through scanpy")
                sc.tl.pca(self.adata, n_comps=self.k)
            self.V = self.adata.varm["PCs"].T
            ev = self.adata.uns["pca"]["variance"]
            self.s = np.sqrt(ev * (self.adata.shape[0] - 1))
            self.U = self.adata.obsm['X_pca'] / self.s[np.newaxis, :]
            # TODO: alternatively use X_pca from harmony (instead of U*s)

    def _setup_covariate_PC_comparison(self):
        # Prep for covariate selection:
        self._calculate_covariate_svd_correlation()
        self._assign_covariates_to_PCs(invert=False)
        self._calculate_covariate_lengths()  # For plotting, later

    def _select_PCs(self, filter_covariate=None, max_k=None):
        # Select SVD components to keep in the correlation calculation
        # TODO: allow invert?
        # TODO: allow multiple components (e.g. celltype + region)
        if filter_covariate is not None and \
                filter_covariate in self.cvlist:
            # Update if it wasn't run originally
            if filter_covariate not in self.cv_pc_ind.keys():
                self._setup_covariate_PC_comparison()
            indices = self.cv_pc_ind[filter_covariate]
            if not self.keep_first_PC:
                indices = indices[indices != 0]
        else:
            indices = None if self.keep_first_PC else \
                np.arange(1, self.U.shape[1])
        if max_k is not None:
            indices = indices[indices <= max_k]
        logging.debug(indices) # TODO: pretty print, maybe #k
        return(indices)

    def _calculate_covariate_svd_correlation(self):
        """Calculate the correlation of covariates with PCs."""
        self.cv_mats = calculate_svd_covar_corr(
            self.U.T, self.adata.obs, self.cvlist, cv_mats=self.cv_mats)

    def _assign_covariates_to_PCs(self, invert=False):
        """Calculate which components are correlated with each covariate."""
        # self._calculate_covariate_svd_correlation()
        for covar in self.cvlist:
            cmat = self.cv_mats[covar].T
            # TODO: fix this normalization issue
            sfact = np.max(self.s) / self.s
            cmat = cmat * sfact[np.newaxis, :]
            cmx = np.max(np.abs(cmat), axis=0)
            # Select indices for each covariate
            if invert:
                # Ones not strongly associated with the covariate
                ind = np.where(cmx < self.cv_cutoff)[0]
            else:
                # Ones strongly associated with the covariate
                ind = np.where(cmx >= self.cv_cutoff)[0]
            self.cv_pc_ind[covar] = ind

    # TODO: Add PC adjustment if necessary
    def _adjust_PCs(self):
        pass

    # TODO: Auxiliary, try to move out of this object
    def _calculate_covariate_lengths(self):
        """Process covariates, for plotting avg heatmap and vs SVD."""
        if not hasattr(self, 'cv_lengths'):
            self.cv_lengths = {}
            self.cv_ticks = {}
        for covar in self.cvlist:
            if covar not in self.cv_lengths.keys():
                cvcol = self.adata.obs[covar]
                if is_categorical_dtype(cvcol):
                    covar_dummy = pd.get_dummies(cvcol)  # To match cv_mats
                    self.cv_lengths[covar] = covar_dummy.shape[1]
                    self.cv_ticks[covar] = covar_dummy.columns.tolist()
                else:
                    self.cv_lengths[covar] = 1
                    self.cv_ticks[covar] = False

    # TODO: put the majority of the correlation function into an import
    def _calculate_correlation(self, indices=None, center=False, power=0,
                               raw=False):
        """
        Calculate the correlation between the genes.
        ---
        Calculates either raw or estimated, using the SVD and selected PCs.
        """
        adjacency = None
        if raw:
            # Raw correlation, centered or not:
            corr = calculate_correlation(self.adata.X, center=center)
        else:
            # Estimate correlation, with a subset of PCs:
            corr = calculate_correlation_estimate(
                self.U, self.s, self.V, power=power, indices=indices)
        return(corr, adjacency)

    def _calculate_correlation_estimate_sd(self, indices=None, power=0):
        # Estimate the std. deviation of the transformed correlation estimate
        # TODO: Remove long-term and replace with bootstraps over batches

        # TODO: Add indices and power variables!
        _, corr_sd = calculate_correlation_estimate_sd(
            self.U, self.s, self.V, nperm=50, seed=self.seed)
        return corr_sd


    # Part II. Functions for graph object, including correlation calculation:
    # -----------------------------------------------------------------------
    # Building adjacency > graph > modules
    # Make graph object > construct_graph does all of these:
    # 5. Build adjacency matrix
    # 6. Filter k-NN (thresholding) > adjacency
    # 6a. Make graph
    # 7. Perform community detection (NMF, BigClam, Leiden)
    # 8. Layout graph
    # 9. Downstream tasks
    # 10. Benchmarking
    def make_graph(self, graph_id, multigraph=False,
                   # For PC selection + correlation estimate:
                   power=0, **kwargs):
        """
        Make graph
        ---
        multigraph: multigraph or single.
        power:      power for eigenvalues
        method:     thresholding method (cutoff, bivariate, sd)
        """
        if multigraph:
            self._make_merged_graph(graph_id, power_list=power, **kwargs)
        else:
            self._make_single_graph(graph_id, power=power, **kwargs)

    # Build the graph object and store in dictionary:
    # TODO: Simplify inheritance of kwargs params:
    def _make_single_graph(self, graph_id,
                           # Correlation options:
                           power=0, filter_covariate=None, raw=False,
                           # Graph options:
                           method='bivariate', # Method for threshold
                           resolution=2,
                           # Processing options:
                           adjacency_only=False,
                           full_graph_only=False,
                           keep_all_z=False,
                           layout=True,
                           **kwargs):
        # 3. Select which PCs to keep:
        indices = self._select_PCs(filter_covariate=filter_covariate)

        # 4. Estimate correlation and 5. Threshold correlation:
        corr, adj = self._construct_adjacency_matrix(
            # Options for constructing correlation
            indices=indices, power=power, raw=raw,
            # Options for thresholding methods:
            method=method,
            **kwargs) # TODO: how to separate out kwargs for this vs. below

        # Make graph with adjacency instead of correlation:
        self._make_single_graph_object(graph_id, corr=corr, adj=adj, **kwargs)

        # Process graph:
        if not adjacency_only and full_graph_only:
            # Compute the full, unaltered, graph for multiplexing,
            # but do not cluster or create modules
            self.graphs[graph_id].construct_graph(
                full_graph=True, modules=False, layout=False)
        else:
            # Build the graph, get modules, and annotate genes:
            # TODO: simplify where adata is called (for populate modules):
            self.graphs[graph_id].construct_graph(
                resolution=resolution, layout=layout)
            self.graphs[graph_id].populate_modules(self.adata, attr='leiden')
            self.graphs[graph_id].match_genes_to_modules(attr='leiden')

    # TODO: Make corr and threshold corr as two different parts:
    def _construct_adjacency_matrix(self,
                                    # Options for correlation:
                                    indices=None, power=0, raw=False,
                                    # Options for thresholding method:
                                    method='bivariate',
                                    # TODO: Add adjacency options here:
                                    **kwargs):
        """
        Construct the adjacency matrix for the given parameters.
        ---
        - Calculates correlation (raw, base)
        - Thresholds the correlation (cutoff, bivariate, sd)
        """
        # 4. Calculate the correlation:
        corr, adjacency = self._calculate_correlation(
            indices=indices, power=power, raw=raw)
        if method == 'sd' and not raw:
            # TODO: update function to use correct indices and power!
            corr_sd = self._calculate_correlation_estimate_sd(
                indices=indices, power=power)
        else:
            corr_sd = None

        # 5. Threshold the correlation using adjacency object:
        adj = adjacency_matrix(
            corr=corr, adjacency=adjacency,
            corr_sd=corr_sd, method=method,
            labels=self.genes, margin=self.margin,
            **kwargs)

        # Return the adjacency object to construct graph
        return(corr, adj)


    # TODO: Simplify inheritance of kwargs params:
    def _make_single_graph_object(self, graph_id, corr,
                                  adj=None, graph=None, # Either one needed
                                  edge_weight=None, min_size=4,
                                  layout_method="fr", **kwargs):
        self.graphs[graph_id] = gene_graph(corr, adj=adj,
                                           graph=graph,
                                           genes=self.genes,
                                           edge_weight=edge_weight,
                                           layout_method=layout_method,
                                           min_size=min_size)

    # Make a merged graph from multiple decorrelation resolutions:
    # TODO: Assign modules later should be using z-scores
    # TODO: If given new z-score, re-threshold all graphs
    def _make_merged_graph(self, graph_id,
                           power_list=[0, .25, .5, .75, 1],
                           filter_covariate=None,
                           resolution=2,
                           keep_all_z=False, # TODO: default to True?
                           **kwargs):
        """Make a merged graph from a list of many parameters."""
        # Make the full list of graphs at each power:
        # TODO: Move make_graphlist into modules?
        # TODO fix, with changes to graphs
        graphlist, graphs = make_graphlist(
            self, power_list=power_list, keep_all_z=keep_all_z,
            filter_covariate=filter_covariate, **kwargs)
        # Multiplex cluster the graphs:
        # NOTE: Leiden partition here, could move to graph utils, add methods
        membership = partition_graphlist(graphlist, resolution=resolution)
        # Calculate the average correlation:
        corr = self._get_average_correlation(graphs)
        # Make the graph object from the union of these graphs:
        # layout and calculate modules
        self._make_merged_graph_object(
            graph_id, graphlist, corr=corr, membership=membership, **kwargs)

    def _get_average_correlation(self, graphs):
        # TODO: Corr is correlation if not loaded as zcutoff.....!
        # NOTE: Extended assignment is always better with zcutoff
        # TODO: Standardize graph: use cutoff to put matrix -> zmat directly
        corr = self.graphs[graphs[0]].corr
        for graph in graphs[1:]:
            # Scale graphs appropriately:
            maxval = np.max(np.abs(self.graphs[graph].corr))
            corr += (self.graphs[graph].corr / maxval)
        corr = corr / (len(graphs) * 1.0)
        return(corr)

    # TODO: should we have this only be the graph object (like single graph
    # object) or full processing, as it is now. If so, rename
    def _make_merged_graph_object(self, graph_id, graphlist, corr,
                                  membership, **kwargs):
        """Combine all graphs, partition modules, layout, and populate."""
        # 1. Combine all of the graphs together:
        graph = self._combine_graphlist(graph_id, graphlist, corr, **kwargs)
        # 2. Partition to modules
        self.graphs[graph_id]._partition_multigraph(
            graph, graphlist, membership, method='leiden')
        # self._partition_multigraph(graph_id, graph, graphlist, membership)
        # 3. Layout the merged graph:
        self.graphs[graph_id].layout_graph()
        # 4. Populate modules with all genes:
        self.graphs[graph_id].populate_modules(self.adata, attr='leiden')
        self.graphs[graph_id].match_genes_to_modules(attr='leiden')

    def _combine_graphlist(self, graph_id, graphlist, corr, **kwargs):
        """Combine all of the graphs in a list of graphs."""
        import igraph
        import contextlib
        # NOTE: Not good! but otherwise igraph prints a huge message.
        with contextlib.redirect_stdout(None):
            # Union and simplify to merge edges
            graph = igraph.union(graphlist)
            graph = graph.simplify(combine_edges="max")

        # Construct the internal graph object for the merged graph:
        self._make_single_graph_object(
            graph_id, corr=corr, graph=graph, **kwargs)
        self.graphs[graph_id].kept_genes = graph.vs['name']
        return(graph)

    # TODO: could update to work with multigraph, or any set of graph options
    def get_k_stats(self, k_list, power=0, resolution=None,
                    raw=False, method='bivariate', filter_covariate=None,
                    **kwargs):
        """Get statistics on # genes and # modules for each k setting."""
        ngenes = []
        nmodules = []
        for k in k_list:
            indices = self._select_PCs(
                filter_covariate=filter_covariate, max_k=k)
            k_graph_id = 'test_k' + str(k)
            # Estimate the correlation and build the adjacency matrix
            corr_subset, adj_subset = self._construct_adjacency_matrix(
                indices=indices, power=power, raw=raw,
                method=method, **kwargs)
            try:
                # Build the graph, - don't layout, but get modules:
                self._make_single_graph_object(
                    k_graph_id, corr=corr_subset, adj=adj_subset, **kwargs)
                self.graphs[k_graph_id].construct_graph(
                    modules=True, resolution=resolution, layout=False)
                # Store number of genes and modules:
                nmodules.append(np.max(
                    self.graphs[k_graph_id].assign['leiden'] + 1))
                ngenes.append(len(self.graphs[k_graph_id].graph.vs))
                # Delete graph after we have metrics:
                del(self.graphs[k_graph_id])
                gc.collect()
            except BaseException:
                ngenes.append(0)
                nmodules.append(0)
            logging.info(f"{k}: ngenes={ngenes} and nmodules={nmodules}")
        return(ngenes, nmodules)


    # Utilities for working directly on a specific graph, maybe not necessary
    # -----------------------------------------------------------------------
    def recluster_graph(self, graph_id, resolution=None):
        """Re-cluster a graph with a different resolution."""
        # NOTE could update to store multiple modules at different resolution
        self.graphs[graph_id].calculate_gene_modules(resolution=resolution)

    def get_modules(self, graph_id, attr="leiden", print_modules=False):
        """Get list of modules from graph and clustering."""
        modules = self.graphs[graph_id].get_modules(
            attr=attr, adata=self.adata, print_modules=print_modules)
        return modules

    def get_module_assignment(self, graph_id, attr="leiden"):
        """Get module assignment for each gene as a pandas DataFrame."""
        mdf = self.graphs[graph_id].get_module_assignment(
            attr=attr, adata=self.adata)
        return mdf

    def find_gene(self, graph_id, gene, return_genes=True, print_genes=True):
        """Find module corresponding to a gene."""
        out = self.graphs[graph_id].find_gene(gene, print_genes=print_genes)
        if return_genes:
            return out

    # Function for saving modules from this object level
    def save_modules(self, graph_id, attr="leiden",
                    as_df=True, filedir="./", filename=None):
        """Save module list for a specific graph as txt or tsv."""
        if filename is None:
            filename = filedir + "module_df_" + self.suffix + "_" + \
                attr + "_" + graph_id
            filename += ".tsv" if as_df else ".txt"

        self.graphs[graph_id].save_modules(
            attr=attr, as_df=as_df, filename=filename)

def calculate_margin_genes(X):
    margin = np.mean(X > 0, axis=0).copy()
    if sparse.issparse(X):
        margin = np.array(margin)[0]
    return(margin)

# NOTE: Had pickling __setstate__ and __getstate__ but removed
