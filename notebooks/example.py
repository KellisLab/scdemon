#!usr/bin/python
"""Example code for using scmodule part of the library."""
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

# Set arguments for plotting outputs:
imgdir = "./"
# Make the scmodule object:
max_k = 100
mod = sm.modules_core(adata,
                      # Tagline for plots
                      suffix=tag,
                      # imgdir=imgdir, TODO: Set global imgdir?
                      # Options for graph creation:
                      estimate_sd=False,
                      svd_k=max_k, filter_expr=0.05, z=4.5,
                      # Overall usage/correlation computation options:
                      calc_raw=False)
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


# Make a graph from multiple resolutions:
# ---------------------------------------
graph_id = 'merge'
mod.make_graph(graph_id, multigraph=True, power=[0,.5,1])

# Plot these:
pl.plot_genes(mod, graph_id, attr="leiden", show_labels=True, width=16)
pl.plot_heatmap_avgexpr(mod, graph_id, cvlist=['leiden'], attr="leiden")
pl.plot_umap_grid(mod, graph_id)


# Make a graph from only specific non-covariate-correlated PCs
# ------------------------------------------------------------
# TODO: cv_cutoff was 2.0
graph_id = 'subset'
# Remove PCs correlated with the cell clustering (leiden)
mod.make_graph(graph_id, filter_covariate="leiden")

# Plot these:
pl.plot_genes(mod, graph_id, attr="leiden", show_labels=True, width=16)
pl.plot_heatmap_avgexpr(mod, graph_id, cvlist=['leiden'], attr="leiden")
pl.plot_umap_grid(mod, graph_id)
