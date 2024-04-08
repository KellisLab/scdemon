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
import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype
import scanpy as sc  # TODO: Don't import if not using adata formulation
from scipy import sparse
import igraph
import contextlib


class modules(object):
    """
    Calculate gene modules from anndata or stand-alone matrix.
    """

    def __init__(self,
                 adata,
                 csuff,
                 h5ad_file=None,  # File corresponding to adata (for saving)
                 imgdir=None,
                 estimate_sd=False,
                 seed=1,
                 filter_expr=None,
                 z=4, svd_k=100,
                 calc_raw=False,
                 use_fbpca=True):
        """Initialize gene modules object."""
        # TODO: implement stand-alone matrix
        # Arguments:
        self.adata = adata
        self.h5ad_file = h5ad_file
        self.genes = self.adata.var_names  # labels
        self.csuff = csuff
        self.z = z
        self.svd_k = svd_k
        self.filter_expr = filter_expr
        self.use_fbpca = use_fbpca
        self.estimate_sd = estimate_sd
        self.seed = seed
        self.calc_raw = calc_raw
        # Directories:
        if imgdir is not None:
            self.imgdir = imgdir
        else:
            self.imgdir = "./"
        # Additional metadata / files:
        self.graphs = {}  # Store computed graphs
        self.cv_mats = {}

    # TODO: Enable this to work with / without an adata object
    def setup(self):
        """Set up standard graph."""
        self._filter_data()
        self._calculate_margin()
        self._calculate_correlation()

    def _filter_data(self):
        # Filter to top expr:
        if self.filter_expr is not None:
            margin = np.mean(self.adata.X > 0, axis=0)
            if sparse.issparse(self.adata.X):
                margin = np.array(margin)[0]
            filtgenes = self.adata.var_names[margin >= self.filter_expr]
            self.adata = self.adata[:, filtgenes]
            self.genes = self.adata.var_names
            print(self.adata.shape)

    def _calculate_margin(self):
        # Calculate margin for filtering later.
        self.margin = np.mean(self.adata.X > 0, axis=0).copy()
        if sparse.issparse(self.adata.X):
            self.margin = np.array(self.margin)[0]

    def _calculate_correlation(self):
        np.random.seed(self.seed)
        if 'X_pca' in self.adata.obsm:
            U = self.adata.obsm["X_pca"]  # TODO: FIX U for scanpy
            s = self.adata.uns["pca"]["variance"]
            V = self.adata.varm["PCs"].T
            self.cobj = correlation(X=self.adata.X,
                                    margin=self.margin,
                                    U=U, s=s, V=V,
                                    k=self.svd_k,
                                    calc_raw=self.calc_raw)

        elif self.use_fbpca:
            # Allow handler to compute PCA with FBPCA library by not feeding in
            # any matrices for U, s, or V
            # TODO: Allow centering if under a certain size / not sparse
            self.cobj = correlation(X=self.adata.X,
                                    margin=self.margin,
                                    k=self.svd_k,
                                    calc_raw=self.calc_raw,
                                    center=False)
        else:
            # Compute PCA with scanpy options if not using FBPCA
            if "X_pca" not in self.adata.obsm.keys():
                logging.debug("Computing PCA through scanpy")
                sc.tl.pca(self.adata, n_comps=self.k)

            U = self.adata.obsm["X_pca"]  # TODO: FIX U for scanpy
            s = self.adata.uns["pca"]["variance"]
            V = self.adata.varm["PCs"].T
            self.cobj = correlation(X=self.adata.X,
                                    margin=self.margin,
                                    U=U, s=s, V=V,
                                    k=self.svd_k,
                                    calc_raw=self.calc_raw)

        # Setup and get the correlation objects
        self.cobj.setup()  # Primarily to calculate PCA if not done yet
        self.corr_est, self.corr_raw = self.cobj.get_correlation()
        # Estimate the std. deviation of the transformed correlation estimate
        if self.estimate_sd:
            self.corr_mean, self.corr_sd = self.cobj.estimate_corr_sd(nperm=50)
        else:
            self.corr_mean, self.corr_sd = None, None

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

    def make_subset_graph(self,
                          graph_id,
                          covar,
                          cv_cutoff=0.4,
                          invert=False,
                          **kwargs):
        """Make a graph from a covariate-correlated subset of SVD dims."""
        if not hasattr(self, "cv_mats") or covar not in self.cv_mats.keys():
            self.calc_svd_corr([covar])
        # TODO: allow multiple components (e.g. celltype + region)
        cmat = self.cv_mats[covar].T
        sfact = np.max(self.cobj.s) / self.cobj.s
        cmat = cmat * sfact[np.newaxis, :]
        cmx = np.max(np.abs(cmat), axis=0)
        if invert:
            k_ind = np.where(cmx < cv_cutoff)[0]
        else:
            k_ind = np.where(cmx >= cv_cutoff)[0]
        nk = len(k_ind)
        logging.info(f"Subsetting to {nk} components for covariate: {covar}")
        self.make_svd_subset_graph(graph_id, k_ind=k_ind, **kwargs)

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

    def _calculate_covariate_lengths(self):
        if not hasattr(self, 'cv_lengths'):
            self.cv_lengths = {}
            self.cv_ticks = {}
        cvlist = self.adata.obs.columns
        for covar in cvlist:
            if covar not in self.cv_lengths.keys():
                cvcol = self.adata.obs[covar]
                if is_categorical_dtype(cvcol):
                    covar_dummy = pd.get_dummies(cvcol)  # To match cv_mats
                    self.cv_lengths[covar] = covar_dummy.shape[1]
                    self.cv_ticks[covar] = covar_dummy.columns.tolist()
                else:
                    self.cv_lengths[covar] = 1
                    self.cv_ticks[covar] = False

    def calc_svd_corr(self, cvlist):
        self.cv_mats = calculate_svd_covar_corr(
            self.cobj.U.T, self.adata.obs, cvlist, cv_mats=self.cv_mats)


    # Functions for saving modules or the full object:
    # ------------------------------------------------
    def save_modules(self, graph_id, attr="leiden", as_df=True,
                     filedir="./", filename=None):
        """Save module list for a specific graph as txt or tsv."""
        if filename is None:
            filename = filedir + "module_df_" + self.csuff + "_" + \
                attr + "_" + graph_id
            filename += ".tsv" if as_df else ".txt"
        self.graphs[graph_id].save_modules(
            attr=attr, as_df=as_df, filename=filename)

    def save_adata(self):
        """Save adata for saving object."""
        # TODO: Currently only saves when file doesn't exist
        # come up with a better choice of when to save adata:
        if not os.path.exists(self.h5ad_file):
            self.adata.write(self.h5ad_file)

    # Functions for properly pickling these objects:
    def __getstate__(self):
        """Get state for object, removing adata, for pickling."""
        # Save the adata object first
        self.save_adata()
        # Copy the object's state from self.__dict__
        state = self.__dict__.copy()
        # Remove the unpicklable entries (adata may be huge):
        del state['adata']
        # TODO: Should we also remove cobj (U, s, V?)
        return state

    def __setstate__(self, state):
        """Load object properly, with adata separately."""
        # Restore instance attributes (i.e., filename and lineno).
        self.__dict__.update(state)
        # Restore the adata with the h5ad_file handle:
        self.adata = sc.read(self.h5ad_file)
