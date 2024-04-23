"""single-cell decorrelated module networks"""
from .modules import modules
from .graph import gene_graph, adjacency_matrix
from .plotting import plot_genes, plot_umap_grid, plot_svd_corr
from .data import snap_colors
from .utils import (
    recipe_preprocess, recipe_annotate, recipe_full,
    get_goterms, vcorrcoef
)

