#!/usr/bin/env python3
"""Gene graph handling."""
# --------------------------------
# Class for handling a gene graph:
# Updated: 06/28/21
# --------------------------------
from ..utils import vcorrcoef
from ..data import snap_colors
from .utils_community_detection import (
    calculate_gene_modules_leiden,
    get_modules_from_partition,
    set_module_colors
)

from .utils_pruning import prune_graph_components, prune_graph_degree

import logging
import numpy as np
import pandas as pd
from scipy import sparse

# For graphs
import igraph


class gene_graph(object):
    def __init__(self,
                 corr,
                 genes,
                 adj=None,
                 graph=None,
                 edge_weight=None,
                 layout_method="fr",
                 min_size=4):
        self.corr = corr
        self.genes = genes
        self.adj = adj      # Adjacency matrix object
        self.graph = graph  # For giving pre-computed graph
        if (graph is None) and (adj is None):
            raise ValueError("Require either graph or adjacency object.")

        # Setup the adjacency matrix processing object:
        self.min_size = min_size
        self.edge_weight = edge_weight
        self.layout_method = layout_method
        self.umat = None

        # Load graph colors resources:
        self.snapcols = snap_colors()
        self.snapcols = list(self.snapcols.to_numpy().T[0])
        self.ncols = len(self.snapcols)

        # For covariate annotations
        self.assign = {}
        self.colors = {}
        # For tracking modules:
        self.modules = {}
        self.mnames = {}
        self.scores = {}
        self.module_match = {}

    # TODO: Add options for different methods of module detection here:
    def construct_graph(self, resolution=2, method='leiden',
                        full_graph=False, modules=True, layout=True):
        """Construct the graph object and find gene modules."""
        # Building adjacency > graph > modules, handled by graph object:
        if self.graph is None:
            # 5. Build adjacency matrix
            self.adjacency, self.kept_genes, self.ind = \
                self.adj.get_adjacency()
            # 6. Create k-NN (thresholding) from adjacency
            # Note: full adjacency graph is for multi-graph purposes.
            self._make_graph(full=full_graph)
        # TODO: check if modules is None, also remake if made new graph
        if modules:
            # 7. Perform community detection (NMF, BigClam, Leiden)
            self.calculate_gene_modules(
                resolution=resolution, method=method)
        if layout:
            # 7a. Layout graph
            self.layout_graph(layout_method=self.layout_method)

    # TODO: Decide if arguments feed in or keep everything in self.
    def _make_graph(self, full=False):
        """Make graph object from adjacency."""
        logging.info("Making graph object from adjacency.")
        if full:
            # Use the full adjacency before adjacency pruning:
            self.adjacency = self.adj.full_adjacency
            self.kept_genes = self.adj.full_labels

        # Get edge list from adjacency matrix:
        self.adjacency = self.adjacency.tocoo()
        edge_list = list(tuple(zip(self.adjacency.row, self.adjacency.col)))
        edge_wt_data = self.adjacency.data

        # Make graph from edge list:
        self.graph = igraph.Graph(edge_list, directed=False,
                                  vertex_attrs={'name': self.kept_genes})
        self.graph.simplify(combine_edges="max")  # Merge edges
        if self.edge_weight is not None:
            self.graph.es["weight"] = self.edge_weight
        else:
            self.graph.es["weight"] = edge_wt_data

        if not full:
            # Clean up very small components (1-2 nodes):
            self.graph = prune_graph_components(
                self.graph, min_size=self.min_size)
            # Clean up orphan nodes with degree 0:
            self.graph = prune_graph_degree(self.graph, degree=0)
            logging.info("Pruned to " + str(len(self.graph.vs)) +
                        " nodes, out of " + str(self.corr.shape[0]))

    def layout_graph(self, layout_method='fr'):
        """Compute graph layout."""
        logging.info("Laying out graph, method=" + str(layout_method))
        if layout_method == "fr":
            self.layout = self.graph.layout_fruchterman_reingold(
                niter=500, grid=False)
        else:
            self.layout = self.graph.layout(layout_method)

    def calculate_gene_modules(self, method="leiden", **kwargs):
        """Calculate modules from gene-gene graph using graph clustering."""
        logging.info("Running module detection using method=" + str(method))
        # TODO: fix assigning method for factors
        if method == 'leiden':
            self.assign[method] = calculate_gene_modules_leiden(
                self.graph, **kwargs)
        else:
            # TODO: implement louvain, resolution seems more tunable?
            raise ValueError("Method '%s' not found" % method)

        # TODO: fix color for factors, maybe just store one color per module
        self.colors[method] = set_module_colors(
            self.assign[method], self.snapcols)

    # TODO: put full handling of multigraph in gene_graph, including graphlist
    def _partition_multigraph(self, graph, graphlist, membership,
                              method = 'leiden'):
        # Turn multiplex partition into modules
        # NOTE: Must reorder due to different order in the merge:
        gn = np.array(graphlist[0].vs['name'])
        reord = np.array([np.where(gn == x)[0][0] for x in graph.vs['name']])
        # gn[reord] == graph.vs['name']  # If we want to check correct.
        membership = np.array(membership)[reord]
        ptns = np.unique(membership)
        partition = [np.where(membership == x)[0] for x in ptns]
        self.assign[method] = get_modules_from_partition(graph, partition)
        self.colors[method] = set_module_colors(
            self.assign[method], self.snapcols)

    # TODO: Allow compute on adjusted adjacency
    def compute_umap(self):
        """Calculate UMAP for correlation estimate underlying graph."""
        import umap
        logging.info("Calculating UMAP for graph adjacency.")
        graph_names = self.graph.vs["name"]
        if self.umat is None or self.umat.shape[1] != len(graph_names):
            # Subset to the appropriate genes for the graph:
            gind = np.array([np.where(self.genes == x)[0][0]
                             for x in graph_names])
            corr_gene_subset = self.corr[gind[:, np.newaxis], gind]
        # Fit the UMAP on the subsetted correlation:
        uw = umap.UMAP()
        self.umat = uw.fit_transform(corr_gene_subset)

    def get_modules(self, attr="leiden", adata=None, print_modules=False):
        """Get list of modules from graph and clustering."""
        if attr not in self.modules.keys():
            # Construct object if necessary:
            self.populate_modules(adata, attr)
        modules = self.modules[attr]
        if print_modules:
            for ll in modules.keys():
                print(ll, " ".join(modules[ll]))
        return modules

    def get_module_assignment(self, attr="leiden", adata=None):
        """Get module assignment for each gene as a pandas DataFrame."""
        if attr not in self.modules.keys():
            # Construct object if necessary:
            self.populate_modules(adata, attr)
        mdf = pd.DataFrame({'gene': self.genes,
                            'module': self.module_match[attr]})
        return mdf

    def find_gene(self, gene, attr='leiden', print_genes=False):
        """Find module corresponding to a gene."""
        mkey = None
        for key in self.modules[attr].keys():
            if gene in self.modules[attr][key]:
                mkey = key
        if mkey is None:
            out = None
            if print_genes:
                print("Gene (" + gene + ") is not present in any modules")
        else:
            out = self.modules[attr][mkey]
            if print_genes:
                print(mkey, " ".join(out))
        return out

    # TODO: Speed this up (possibly avg. by tform multiplication)
    # Modules lists functions (requires external adata or X):
    # TODO: Update to use X only, not adata.X (fix subsetting)
    def populate_modules(self, adata, attr='leiden'):
        """Populate modules data."""
        if adata is None:
            raise TypeError("adata is None, need adata to populate modules")
        logging.info("Populating modules data")
        # TODO: make orderedDict for the module names?
        partition = self.assign[attr]
        nam = np.array(self.graph.vs["name"])
        modules = np.unique(partition)
        issp = sparse.issparse(adata.X)  # TODO: get this into graph class
        # Initialize:
        self.modules[attr] = {}
        self.mnames[attr] = [""] * (np.max(modules) + 1)
        self.scores[attr] = np.zeros((len(adata), np.max(modules) + 1))
        for i in modules:
            mgenes = np.sort(nam[np.where(partition == i)[0]])
            subadata = adata[:, mgenes]
            avgexpr = np.mean(subadata.X, axis=1)
            if issp:
                avgexpr = np.array(avgexpr).T[0]
            else:
                avgexpr = avgexpr.T
            self.modules[attr][i] = mgenes
            self.scores[attr][:, i] = avgexpr  # For plotting on cell umap
            # Find the top 3 genes (for naming module):
            cv = vcorrcoef(subadata.X.T, avgexpr)
            topgenes = ", ".join(mgenes[np.argsort(-cv)][0:3])
            ngstr = str(len(mgenes))
            ngstr = ngstr + " gene" if ngstr == "1" else ngstr + " genes"
            self.mnames[attr][i] = "M" + str(i) + \
                " (" + ngstr + "): " + topgenes

    # TODO: Speed this up - sort step can be slow
    def match_genes_to_modules(self, attr='leiden'):
        """Match genes to the closest to module by its correlation."""
        logging.info("Matching genes to modules")
        # Nearest module by top correlations across the full dataset:
        if attr not in self.module_match.keys():
            self.module_match[attr] = np.zeros(
                (self.corr.shape[0], len(self.modules[attr])))
            mxval = np.max(np.abs(self.corr)) * 10  # Ensure core stays
            for key in self.modules[attr].keys():
                # Subset correlation to all genes vs. module genes:
                genes = self.modules[attr][key]
                lind = np.in1d(self.genes, genes)
                subcorr = self.corr[:, lind]
                # Score a module by the average corr. of top 5 module genes
                self.module_match[attr][:, key] = np.mean(
                    np.sort(subcorr, 1)[:, -5:], 1)
                # Ensure core genes in module stay in it in the full assign:
                self.module_match[attr][lind, key] = mxval
            # Set each gene's module to the top average correlation:
            self.module_match[attr] = np.argmax(self.module_match[attr], 1)

    # Add save functions here (graph and full object as well):
    def save_modules(self, attr="leiden", as_df=True,
                     filename=None, adata=None):
        """Save module list as txt or tsv."""
        if as_df:
            mdf = pd.DataFrame({"gene": np.array(self.graph.vs["name"]),
                                "leiden": self.assign[attr],
                                "color": self.colors[attr]})
            mdf.to_csv(filename, sep="\t")
        else:
            mod_list = self.get_modules(
                attr=attr, adata=adata, print_modules=False)
            with open(filename, "w") as f:
                for key in mod_list.keys():
                    line = str(key) + ": " + " ".join(mod_list[key]) + "\n"
                    f.write(line)
        logging.info("Wrote modules to:" + str(filename))

