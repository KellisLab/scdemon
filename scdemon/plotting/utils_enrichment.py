#!/usr/bin/env python3
"""Functions for module enrichments on dataFrames."""
import numpy as np
import pandas as pd


def calc_df_enrichment(df, col, genes, module_match):
    """
    Calculate enrichment of modules in dataframe.

    Given a dataframe with genes (df.gene) and some split (df[col]),
    calculate the enrichment of each split in the modules
    """
    from scipy.stats import hypergeom
    # Count genes in each:
    keys = np.flip(np.sort(pd.unique(df[col])))
    cmat = np.zeros((np.max(module_match) + 1, len(keys)), int)
    for i, key in enumerate(keys):
        df_genes = df.gene[df[col] == key]
        ind = np.in1d(genes, df_genes)
        u, c = np.unique(module_match[ind], return_counts=True)
        cmat[u, i] = c

    # Hypergeom for each box:
    rmarg = np.sum(cmat, 1)
    cmarg = np.sum(cmat, 0)
    pmat = np.zeros(cmat.shape)
    rmat = np.zeros(cmat.shape)
    M = np.sum(cmat)  # Total # genes
    stats = []
    for i in range(cmat.shape[0]):
        for j in range(cmat.shape[1]):
            x = cmat[i, j]  # Drawn matching module
            n = rmarg[i]  # Any matching module
            N = cmarg[j]  # Total drawn
            pval = hypergeom.sf(x - 1, M, n, N)
            ratio = (x / N) / (n / M)
            pmat[i, j] = -np.log10(pval)
            rmat[i, j] = ratio
            # Track statistics to return as DataFrame:
            stats.append([i, j, x, n, N, M, ratio, pval])

    # Statistics table:
    statsdf = pd.DataFrame(
        np.array(stats),
        columns=[
            'module',
            'j',
            'x',
            'n',
            'N',
            'M',
            'ratio',
            'p'])
    for column in ['module', 'j', 'x', 'n', 'N', 'M']:
        statsdf[column] = statsdf[column].astype(int)

    statsdf['key'] = keys[statsdf.j.to_numpy()]

    return pmat, rmat, keys, statsdf


def format_enrichment_matrices(rmat, pmat):
    """Format the enrichment matrices for plotting."""
    # Set NaNs and 0 to min (for plotting)
    rmat[np.isnan(rmat)] = 1
    rmat[rmat == 0] = np.min(rmat[rmat != 0])
    # Format the annotation for plotting:
    annot = np.empty(rmat.shape, dtype="<U10")
    annot[:] = ""
    annot[10 ** -pmat < 0.05] = "*"
    annot[10 ** -pmat < 0.01] = "**"
    annot[10 ** -pmat < 0.001] = "***"
    return (rmat, annot)

