#!/usr/bin/env python3
"""Plotting genes on the UMAP basis or the gene-gene graph."""
import os
import re
import gc
import logging
import numpy as np
import pandas as pd
import scanpy as sc  # Need for genes' logFC vs. a covariate

# Plotting:
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize, rgb2hex
import seaborn as sns

from ..graph import gene_graph
from .utils_graph_plotting import (
    _plot_gene_graph_core, _plot_gene_umap_core
)



# Wrapper for plotting either from overall object or from graph:
def plot_genes(obj, graph_id=None, basis='graph', attr=None,
               color=None, imgdir='./', suffix=None, ext='png', **kwargs):
    """\
        Plot all genes in a graph on a graph or UMAP basis

        Parameters
        ----------
        obj : gene_graph | modules
            Object (``gene_graph`` or ``modules``) with modules to plot
        graph_id : str
            Name of graph to work with
        attr : str
            Name of modules in the graph
        basis : str
            Plotting basis, either ``'graph'`` or ``'umap'``
        color : str | list
            Color(s) override
        imgdir : str
            Directory for images
        suffix : str
            Suffix for image filenames
        ext : str
            Extension for image (``'png'`` or ``'pdf'``)
        **kwargs
            Additional arguments for ``_plot_genes_for_graph``, including ``width``, ``title``, and ``ax``
    """
    # Check if object is a graph or not:
    if type(obj) is gene_graph:
        _plot_genes_for_graph(obj, basis=basis, attr=attr, color=color,
                              imgdir=imgdir, suffix=suffix, ext=ext, **kwargs)
    else:
        # If not a graph, check for the graph object underneath it:
        if (graph_id is None) or (not hasattr(obj, 'graphs')) or \
                (graph_id not in obj.graphs):
            raise ValueError(
                "Graph with ID '%s' doesn't exist in this object." % graph_id)
        else:
            # If the graph exists, plot the genes under this graph:
            # NOTE: Also plots it to the object's imgdir
            if suffix is None:
                suffix = obj.suffix + "_" + graph_id
            _plot_genes_for_graph(
                obj.graphs[graph_id], basis=basis, attr=attr, color=color,
                imgdir=imgdir, suffix=suffix, ext=ext, **kwargs)


def _plot_genes_for_graph(obj, basis='graph', attr=None, color=None,
                          imgdir='./', suffix=None, ext='png', **kwargs):
    """
    Plot the genes under a specific basis (umap/graph).

    obj: Object is a graph, usually from modules, mod_obj.graphs[graph_id]
    """
    # Set color from attribute and specify plotname:
    plotname = imgdir + basis + "_"
    if attr is None:
        col = None
    else:
        if color is not None:
            obj.colors[attr] = color
        col = obj.colors[attr]
        plotname += str(attr) + "_"
    if suffix is not None:
        plotname += suffix + "." + ext

    # Plot the genes on a specific basis (graph/umap):
    if basis == 'graph':
        _plot_gene_graph_from_gene_graph(
            obj, col=col, plotname=plotname, attr=attr, **kwargs)
    elif basis == 'umap':
        if obj.umat is None:
            obj.compute_umap()
        _plot_gene_umap_from_gene_graph(
            obj, col=col, plotname=plotname, **kwargs)


# Calculate UMAP representation and plot the genes on it:
def _plot_gene_umap_from_gene_graph(obj, col=None, plotname=None,
                                    ax=None, title=None, width=12):
    """Plot the genes on the UMAP basis with a specified attribute."""
    _plot_gene_umap_core(obj.umat, col=col,
                         width=width, title=title,
                         plotname=plotname, ax=ax)


def _plot_gene_graph_from_gene_graph(obj, col=None, width=24, plotname=None,
                                     ax=None, title=None, attr='leiden',
                                     show_labels=True, frac_labels=1.0,
                                     adjust_labels=False):
    """Plot the genes on the gene-gene graph with a specified attribute."""
    if not hasattr(obj, 'layout'):
        obj.layout_graph(layout_method=obj.layout_method)

    _plot_gene_graph_core(graph=obj.graph, layout=obj.layout,
                          assign=obj.assign['leiden'],
                          plotname=plotname, ax=ax, title=title,
                          col=col, width=width,
                          show_labels=show_labels,
                          frac_labels=frac_labels,
                          adjust_labels=adjust_labels)


# TODO: Improvements here:
# - Add plot of arbitrary vectors (e.g. margin)
# - Figure out to deal with non-categorical
# - Separate out the gene logFC score compute from the plotting
# - Add colorbar
def plot_gene_logfc(obj, graph_id, basis='graph', attr="celltype",
                    fc_cut=2, p_cut=0.05, cmap = plt.cm.RdYlBu,
                    imgdir='./', **kwargs):
    """\
        Plot logFC for all genes against a specific covariate on a graph or UMAP basis

        Parameters
        ----------
        obj : modules
            Object (``modules``) with graph and modules to plot
        graph_id : str
            Name of graph to work with
        attr : str
            Name of modules in the graph
        basis : str
            Plotting basis, either ``'graph'`` or ``'umap'``
        fc_cut : float
            Cap for top and bottom of logFC scale
        p_cut : float
            p-value cutoff for the logFC values
        cmap
            matplotlib colormap to use (default ``plt.cm.RdYlBu``)
        imgdir : str
            Directory for images
        **kwargs
            Additional arguments for ``plot_genes``, including ``width``, ``title``, and ``ax``
    """
    logging.info("Running rank genes on '%s'" % attr)
    sc.tl.rank_genes_groups(obj.adata, attr)
    iscat = obj.adata.obs[attr].dtype.name == "category"
    # Color map:
    norm = Normalize(vmin=-fc_cut, vmax=fc_cut)
    # Turn into assignment:
    narr = obj.adata.uns["rank_genes_groups"]["names"]
    larr = obj.adata.uns["rank_genes_groups"]["logfoldchanges"]
    parr = obj.adata.uns["rank_genes_groups"]["pvals_adj"]
    # Graph object:
    graph_object = obj.graphs[graph_id]
    names = list(narr.dtype.names)
    for name in names:
        namestr = attr + "_" + re.sub("[ /]", "_", name)
        logging.info("Plotting for %s" % namestr)
        n1 = np.array(narr[name])
        xind = np.array([np.where(n1 == x)[0][0]
                            for x in graph_object.graph.vs["name"]])
        n1 = n1[xind]
        l1 = np.array(larr[name])[xind]
        p1 = np.array(parr[name])[xind]
        l1 = l1 * (p1 < p_cut)
        l1[l1 < -fc_cut] = -fc_cut
        l1[l1 > fc_cut] = fc_cut
        # Color mapping:
        vcols = [rgb2hex(c) for c in cmap(norm(-l1))]
        # Plot graph, passing colors to genes:
        plot_genes(graph_object, basis=basis,
                   attr=namestr, color=vcols,
                   suffix=obj.suffix + "_" + graph_id,
                   imgdir=imgdir, **kwargs)
