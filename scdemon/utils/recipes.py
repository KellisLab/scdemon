#!/usr/bin/env python3
"""Simple recipes for scanpy pre-processing."""
import numpy as np
import scanpy as sc
import logging


def recipe_preprocess(adata, min_cells=3, min_genes=100):
    """Simple scanpy recipe for single-cell preprocessing for examples"""
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    logging.info("Preprocessed example dataset")


def recipe_annotate(adata):
    """Simple scanpy recipe for single-cell annotation for examples"""
    sc.tl.pca(adata, n_comps=100)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    logging.info("Annotated example dataset")


def recipe_full(adata, preprocess=True, annotate=True,
                min_cells=3, min_genes=100):
    """Simple scanpy recipe for single-cell preprocessing and annotation for examples"""
    if preprocess:
        recipe_preprocess(adata, min_cells=min_cells, min_genes=min_genes)
    if annotate:
        recipe_annotate(adata)
