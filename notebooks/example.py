#!/usr/bin/env python3
"""Example code for using scdemon part of the library."""
# ---------------------------------------------------------------------
# Example - compute co-expression modules using a scanpy anndata object
# Updated: 04/08/24
# ---------------------------------------------------------------------
import logging
import os
import scanpy as sc
import numpy as np

import scdemon as sm
from scdemon.utils import recipe_full
from scdemon import plotting as pl

# Set logging level:
logging.basicConfig(level=logging.INFO)


# Load and prep datasets:
# -----------------------
# Or one of scanpy's datasets:
tag = "pbmc_example"
annfile = tag + "_ann.h5ad"
if os.path.exists(annfile):
    adata = sc.read(annfile)
else:
    adata = sc.datasets.pbmc3k()
    recipe_full(adata, preprocess=True, annotate=True)
    adata.write(annfile)

logging.info("Loaded example dataset")

# Make the modules handling object:
max_k = 100
mod = sm.modules_core(adata, suffix=tag, # Tagline for plots
                      # Options for graph creation:
                      svd_k=max_k, filter_expr=0.05)
mod.setup()  # Setup the object


# Make graph using the selected parameters for basic analysis:
# ------------------------------------------------------------
graph_id = "base"
mod.make_graph(graph_id, resolution=2.5)

# Plot genes on gene-gene graph and on gene UMAP basis
pl.plot_genes(mod, graph_id, attr="leiden", show_labels=True, width=16)
pl.plot_genes(mod, graph_id, basis='umap', attr="leiden", width=16)

# Plot module expression on the cell UMAP basis:
pl.plot_umap_grid(mod, graph_id)

# Get the modules/print out:
mlist = mod.get_modules(graph_id, print_modules=False)
mod.save_modules(graph_id)

# Get functional enrichments for the modules:
gpres = sm.get_goterms(mod, graph_id)


# We can plot additional plots on the gene modules graph or UMAP:
# ---------------------------------------------------------------
# Plot the logFC for a specific covariate on the graph:
graph_id = "base"
covariate = "leiden"
pl.plot_gene_logfc(mod, graph_id, attr=covariate,
                   show_labels=False, width=16, fc_cut=2)

# Plot correlation of SVD components with covariates:
cvlist = ["n_genes", "leiden"]
pl.plot_svd_corr(mod, cvlist)

# Plot the average expression of leiden modules on covariates:
pl.plot_heatmap_avgexpr(mod, graph_id, cvlist=cvlist, attr="leiden")


# Make some other graphs based on different methods:
# --------------------------------------------------
# 1. Make a graph from multiple resolutions:
mod.make_graph('merge', multigraph=True, power=[0,.5,1])

# 2. Make a graph from the raw correlation:
mod.make_graph('raw', raw=True)

# 3. Remove PCs correlated with the cell clustering (leiden)
mod.make_graph('subset', filter_covariate="leiden")


# Plot these graphs:
# ------------------
graphlist = ['merge', 'subset', 'raw']
for graph_id in graphlist:
    pl.plot_genes(mod, graph_id, attr="leiden", show_labels=True, width=16)
    pl.plot_heatmap_avgexpr(mod, graph_id, cvlist=['leiden'], attr="leiden")
    pl.plot_umap_grid(mod, graph_id)


# TODO: Add get_k_stats test too

