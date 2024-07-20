#!/usr/bin/env python3
"""Example code for using scdemon part of the library."""
# ---------------------------------------------------------------------
# Example - compute co-expression modules using a scanpy anndata object
# Updated: 04/08/24
# ---------------------------------------------------------------------
import os
import logging
import numpy as np
import pandas as pd

import scanpy as sc
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
mod = sm.modules(adata, suffix=tag, # Tagline for plots
                 # Options for graph creation:
                 k=max_k, filter_expr=0.05)
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

# Get the modules as dictionary of lists, and/or print them out:
mlist = mod.get_modules(graph_id, print_modules=True)
# Alternative format: module assignments in DF
moddf = mod.get_module_assignment(graph_id)
# Save them to a file:
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
covariate_list = ["n_genes", "leiden"]
pl.plot_svd_corr(mod, covariate_list)

# Plot the average expression of leiden modules on covariates:
pl.plot_heatmap_avgexpr(mod, graph_id, covariate_list=covariate_list, attr="leiden")


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
    pl.plot_heatmap_avgexpr(mod, graph_id, covariate_list=['leiden'], attr="leiden")
    pl.plot_umap_grid(mod, graph_id)


# Other functionalities:
# ----------------------
# Find the module for a specific gene:
mod.find_gene(graph_id, 'CD74')

# Recluster an existing graph:
mod.recluster_graph(graph_id, resolution=2.5)

# Test module recovery for varying k truncation for SVD:
klist = list(np.arange(30, max_k, 10)) + [max_k]
ng, nm = mod.get_k_stats(k_list=klist, power=0)
score = np.array(ng) / 100 + np.array(nm)
paramdf = pd.DataFrame({'k': klist, 'ng': ng, 'nm': nm, 'score': score})

