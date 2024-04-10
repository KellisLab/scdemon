#!/usr/bin/python
"""[WORKING] Reduced class for computing modules on adata object."""

# Internal:
from .graph import gene_graph
from .correlation import correlation
from .auxiliary import calculate_svd_covar_corr
from .utils_multigraph import make_graphlist, partition_graphlist

import os
import re
import gc
import logging
import contextlib

import fbpca
import igraph
import numpy as np
import pandas as pd
from scipy import sparse
from pandas.api.types import is_categorical_dtype

import scanpy as sc
from anndata import AnnData


# Removed options
# csuff,
# imgdir=None,
# h5ad_file=None,  # File corresponding to adata (for saving)

# TODO: standardize SVD vs. PCs nomenclature

def run_modules(
    adata,
    estimate_sd=False,
    seed=1, # TODO: add fixed random seed on construction
    filter_expr=None,
    z=4,
    svd_k=100,
    calc_raw=False,
    use_fbpca=True,
):
    """Run modules"""
    # TODO: Farm out instances:
    if type(obj) is AnnData:
        # Create modules object differently depending on object.
        mod = modules(adata,
                      csuff=csuff,
                      estimate_sd=estimate_sd,
                      seed=seed,
                      filter_expr=filter_expr,
                      z=z,
                      svd_k=svd_k,
                      calc_raw=calc_raw
                      use_fbpca=use_fbpca)
    else:
        # NOTE: If no anndata, where do pieces go?
        pass
    mod.setup()  # Setup the object
    # TODO: how to handle different graphs?
    graph_id = "base"
    # clustering resolution to main call
    mod.make_graph(graph_id, resolution=2.5)
    # TODO: Enable this to work with / without an adata object
    # Maybe two parts - one is setup > correlation, second is make_graph >  modules


class modules(object):
    """
    Calculate gene modules from anndata or stand-alone matrix.
    """

    def __init__(self,
                 adata,
                 X=None, meta=None, # TODO: Add alternative if no adata.obs
                 U=None, s=None, V=None, # If we want to pass these in

                 seed=1,
                 estimate_sd=False,
                 filter_expr=None,
                 svd_k=100, # TODO: change to k, change example code to match
                 calc_raw=False,
                 use_fbpca=True
                 cv_cutoff=0.4,
                 # Arguments for graph:
                 z=4,
                 ):
        """Initialize gene modules object."""
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

        # Arguments for setup steps (up to correlation estimate)
        self.k = svd_k
        self.filter_expr = filter_expr
        self.use_fbpca = use_fbpca
        self.seed = seed
        self.calc_raw = calc_raw
        self.estimate_sd = estimate_sd
        self.cv_cutoff = cv_cutoff

        # Arguments for graph building steps:
        self.z = z

        # Additional metadata / files:
        self.graphs = {}  # Store computed graphs
        self.cv_mats = {} # For correlation of covariates vs. PCs
        self.cv_pc_ind = {} # For filtering PCs by covariate

        # Initialize with the random seed
        np.random.seed(self.seed)

    # TODO: Allow options:
    # TODO: Possibly separate into setup + calculate correlation,
    # so can we can run in multiple times, make multiple graphs
    def setup(self):
        """Set up the dataset and get the correlation estimate."""
        # 1a. normalize (done outside of this)
        self._filter_by_genes() # 1b. subset to genes above expr cutoff
        self._calculate_PCA() # 2. perform SVD or get pre-computed SVD
        self._select_PCs() # 2a. filter components
        self._adjust_PCs() # 3. PC adjustment # TODO: Put before or after 2a?
        self._calculate_correlation() # 4. Estimate correlation

    def build_graph(self):
        pass

    # Building adjacency > graph > modules
    # 5. Build adjacency matrix
    # 6. Filter k-NN (thresholding)
    # 7. Perform community detection (NMF, BigClam, Leiden)
    # 8. Downstream tasks
    # 9. Benchmarking

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

    def _calculate_PCA(self):
        """
        Calculate PCA components U, s, V if they are not already computed.
        ---
        Defaults to X_pca in adata.obsm if available.
        Otherwise uses fbpca or sc.tl.pca
        """
        if not self._check_PCA():
            if self.use_fbpca and not 'X_pca' in self.adata.obsm:
                # TODO: Check with Ben if other preferred PCA method
                logging.info("Calculating PCA of X with fbpca")
                self.U, self.s, self.V = fbpca.pca(
                    self.X, k=self.k, raw=not self.center)
            else:
                # Compute PCA with scanpy options if not using fbpca
                if "X_pca" not in self.adata.obsm.keys():
                    logging.debug("Computing PCA through scanpy")
                    sc.tl.pca(self.adata, n_comps=self.k)
                self.U = self.adata.obsm["X_pca"]  # TODO: FIX U for scanpy
                self.s = self.adata.uns["pca"]["variance"]
                self.V = self.adata.varm["PCs"].T

    # TODO: Make this optional
    def _select_PCs(self):
        if self.filter_covariate is not None:
            self._calculate_covariate_svd_correlation()
            self._assign_covariates_to_PCs(invert=False)
            self._calculate_covariate_lengths()  # For plotting, later
            # Select SVD components to keep
            # TODO: allow multiple components (e.g. celltype + region)
            covar = self.filter_covariate[covar]
            self.kept_ind = self.cv_pc_ind[covar]

    def _calculate_covariate_svd_correlation(self):
        """Calculate the correlation of covariates with PCs."""
        self.cv_mats = calculate_svd_covar_corr(
            self.U.T, self.adata.obs, self.cvlist, cv_mats=self.cv_mats)

    # TODO: rename.
    def _assign_covariates_to_PCs(self, invert=False):
        """Calculate which components are correlated with each covariate."""
        self.calc_svd_corr([covar])
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

    def _calculate_correlation(self, center=False):
        # Create correlation object:
        # TODO: Change to allow centering if > set size or already dense
        self.cobj = correlation(X=self.adata.X,
                                margin=self.margin,
                                U=U, s=s, V=V,
                                k=self.svd_k,
                                calc_raw=self.calc_raw,
                                center=center)

        # Calculate correlation:

        # TODO: remove getting raw correlation if not necessary
        self.corr_est, self.corr_raw = self.cobj.get_correlation()

        # TODO: remove this
        # Estimate the std. deviation of the transformed correlation estimate
        if self.estimate_sd:
            self.corr_mean, self.corr_sd = self.cobj.estimate_corr_sd(
                nperm=50)
        else:
            self.corr_mean, self.corr_sd = None, None

        # If subset instead:
        if self.filter_covariate is not None:
            # corr_subset = self.cobj.get_correlation_subset(k_ind, power=power)



    # TODO: Simplify inheritance of kwargs params:
    def _make_graph_object(self, graph_id, corr,
                           graph=None, corr_sd=None,
                           cutoff=None, edge_weight=None,
                           knn_k=None, scale=None,
                           z=None, degree_cutoff=0,
                           min_size=4, use_zscore=True,
                           layout_method="fr"):
        if z is not None:
            self.z = z
        self.graphs[graph_id] = gene_graph(
            corr,
            genes=self.genes,
            graph=graph,
            corr_sd=corr_sd,
            cutoff=cutoff,
            use_zscore=use_zscore,
            edge_weight=edge_weight,
            margin=self.margin,
            layout_method=layout_method,
            knn_k=knn_k,
            z=self.z,
            scale=scale,
            min_size=min_size,
            degree_cutoff=degree_cutoff)

    # Make a merged graph from multiple decorrelation resolutions:
    # NOTE: Performs better than each graph singly at calling modules
    # TODO: Figure out how to merge...
    # TODO: Assign modules later should be using z-scores
    # TODO: If given new z-score, re-threshold all graphs
    def make_merged_graph(self,
                          graph_id,
                          power_list=[0, .25, .5, .75, 1],
                          resolution=2,
                          keep_all_z=False,
                          **kwargs):
        """Make a merged graph from a list of many parameters."""
        # Make the full list of graphs at each power:
        # TODO: Move make_graphlist into modules?
        graphlist, graphs = make_graphlist(self, plist=power_list,
                                           keep_all_z=keep_all_z)
        # Multiplex cluster the graphs:
        membership = partition_graphlist(graphlist, resolution=resolution)
        # Calculate the average correlation:
        # TODO: Corr is correlation if not loaded as zcutoff.....!
        # NOTE: Extended assignment is always better with zcutoff
        # TODO: Standardize graph: use cutoff to put matrix -> zmat directly
        corr = self.graphs[graphs[0]].corr
        for graph in graphs[1:]:
            # Scale graphs appropriately:
            maxval = np.max(np.abs(self.graphs[graph].corr))
            corr += (self.graphs[graph].corr / maxval)
        corr = corr / (len(graphs) * 1.0)

        # Make the graph object with the pre-computed graph:
        # NOTE: Not good! but otherwise igraph prints a huge message...
        with contextlib.redirect_stdout(None):
            graph = igraph.union(graphlist)
        graph = graph.simplify(combine_edges="max")  # Merge edges
        self._make_graph_object(
            graph_id, corr=corr, graph=graph, **kwargs)
        self.graphs[graph_id].kept_genes = graph.vs['name']

        # Turn multiplex partition into modules (reorder due to merge):
        gn = np.array(graphlist[0].vs['name'])
        reord = np.array([np.where(gn == x)[0][0] for x in graph.vs['name']])
        # gn[reord] == graph.vs['name']  # If we want to check correct.
        membership = np.array(membership)[reord]
        ptns = np.unique(membership)
        partition = [np.where(membership == x)[0] for x in ptns]
        self.graphs[graph_id].get_modules_from_partition(partition, 'leiden')

        # Layout the merged graph:
        self.graphs[graph_id].layout_graph()
        # Populate modules??
        self.graphs[graph_id].populate_modules(self.adata, attr='leiden')
        self.graphs[graph_id].match_genes_to_modules(attr='leiden')

    # Build the graph object and store in dictionary:
    # TODO: Clarify the hierarchy of options - merge use_zscore and z and
    # cutoff if necessary.
    def make_graph(self,
                   graph_id,
                   corr=None,
                   corr_sd=None,
                   power=0,
                   resolution=2,
                   layout=True,
                   adjacency_only=False,
                   full_graph_only=False,
                   keep_all_z=True,
                   **kwargs):
        # For debugging purposes (testing power, if recomputing, etc.)
        self.corr_est, self.corr_raw = self.cobj.get_correlation(power=power)
        # Set parameters + correlations:
        # TODO: Check kwargs passed properly (corr, corr_sd?)
        if corr is None:
            if graph_id == "raw":
                if self.corr_raw is None:
                    self.corr_raw = self.cobj.get_raw_correlation(force=True)
                usecorr = self.corr_raw
                usesd = None
            else:
                usecorr = self.corr_est
                usesd = self.corr_sd
        else:
            usecorr = corr
            usesd = corr_sd
        # TODO: Simplify inheritance of kwargs params:
        # Make the actual graph object:
        self._make_graph_object(
            graph_id, corr=usecorr, corr_sd=usesd, **kwargs)
        # Process graph:
        if adjacency_only:
            # Only compute adjacency, for bootstraps
            _, _, _ = self.graphs[graph_id].adj.get_adjacency()
        elif full_graph_only:
            # Compute the full, unaltered, graph for multiplexing
            self.graphs[graph_id].construct_full_graph(keep_all_z=keep_all_z)
        else:
            # Build the graph, get modules, and annotate genes:
            # TODO: simplify where adata is called (for populate modules):
            self.graphs[graph_id].construct_graph(
                resolution=resolution, layout=layout)
            self.graphs[graph_id].populate_modules(self.adata, attr='leiden')
            self.graphs[graph_id].match_genes_to_modules(attr='leiden')

    def make_svd_subset_graph(self, graph_id, k_ind, power=0, **kwargs):
        """Make a graph from a subset of SVD dims."""
        # Make corr matrix:
        corr_subset = self.cobj.get_correlation_subset(k_ind, power=power)
        self.make_graph(graph_id, corr=corr_subset, **kwargs)


    def get_k_stats(self, k_list, power=0, resolution=None, **kwargs):
        """Get statistics on # genes and # modules for each k setting."""
        ngenes = []
        nmodules = []
        for k in k_list:
            ind = np.arange(k)
            graph_id = 'test_k' + str(k)
            # Build the graph - don't layout, but get modules:
            corr_subset = self.cobj.get_correlation_subset(ind, power=power)
            try:
                self._make_graph_object(graph_id, corr=corr_subset, **kwargs)
                self.graphs[graph_id].construct_graph(
                    resolution=resolution, layout=False)
                # Store number of genes and modules:
                nmodules.append(np.max(
                    self.graphs[graph_id].assign['leiden'] + 1))
                ngenes.append(len(self.graphs[graph_id].graph.vs))
                # Delete graph after we have metrics:
                del(self.graphs[graph_id])
                gc.collect()
            except BaseException:
                ngenes.append(0)
                nmodules.append(0)
            logging.info(f"{k}: ngenes={ngenes} and nmodules={nmodules}")
        return(ngenes, nmodules)


    # TODO: Put utilities for graph object, somewhere more reasonable
    def recluster_graph(self, graph_id, resolution=None):
        """Re-cluster graph with different resolution."""
        # TODO: Allow multiple modules at different resolution
        self.graphs[graph_id].calculate_gene_modules(resolution=resolution)

    def get_modules(self, graph_id, attr="leiden", print_modules=False):
        """Get list of modules from graph and clustering."""
        modules = self.graphs[graph_id].get_modules(
            attr, print_modules=print_modules)
        return modules

    def get_module_assignment(self, graph_id, attr="leiden"):
        """Get module assignment for each gene as a pandas DataFrame."""
        mdf = self.graphs[graph_id].get_module_assignment(attr=attr)
        return mdf

    def find_gene(self, graph_id, gene, return_genes=True, print_genes=True):
        """Find module corresponding to a gene."""
        out = self.graphs[graph_id].find_gene(gene, print_genes=print_genes)
        if return_genes:
            return out



def calculate_margin_genes(X):
    margin = np.mean(X > 0, axis=0).copy()
    if sparse.issparse(X):
        margin = np.array(margin)[0]
    return(margin)

# Functions for saving modules:
def save_modules(obj, graph_id, suffix, attr="leiden",
                 as_df=True, filedir="./", filename=None):
    """Save module list for a specific graph as txt or tsv."""
    if filename is None:
        filename = filedir + "module_df_" + suffix + "_" + \
            attr + "_" + graph_id
        filename += ".tsv" if as_df else ".txt"

    # TODO: Move this function out of graphs
    obj.graphs[graph_id].save_modules(
        attr=attr, as_df=as_df, filename=filename)

# NOTE: Had pickling __setstate__ and __getstate__ but removed
