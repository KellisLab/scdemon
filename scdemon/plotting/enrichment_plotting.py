#!/usr/bin/env python3
"""Plot the enrichment of covariates on graphs."""
import logging
import numpy as np
import pandas as pd

# Plotting:
from matplotlib import pyplot as plt
import seaborn as sns

from .utils_enrichment import calc_df_enrichment, format_enrichment_matrices


def plot_df_enrichment(obj, graph_id, df, col, suffix, imgdir='./',
                       attr='leiden', title=None, ext="png"):
    """\
        Calculate and plot enrichment of sets of genes (e.g. DEGs) against modules

        Parameters
        ----------
        obj : gene_graph | modules
            Object (``gene_graph`` or ``modules``) with modules to plot
        graph_id : str
            Name of graph to work with
        df : pandas.DataFrame
            Dataframe containing gene sets, with gene column named ``'gene'``
        col : str
            Column name in ``df`` that gives the gene sets
        suffix : str
            Suffix for the plot filename
        attr : str
            Name of modules in the graph
        title : str
            Title for plot
        imgdir : str
            Directory for images
        ext : str
            Extension for file (either ``'pdf'`` or ``'png'``)
    """
    # Check if object is a graph or not:
    if type(obj) is gene_graph:
        graph_obj = obj
        plotname = imgdir + "modules_degenr_heatmap_" + suffix + "." + ext
    else:
        if (not hasattr(obj, 'graphs')) or (graph_id not in obj.graphs):
            raise ValueError(
                "Graph with ID '%s' doesn't exist in this object." % graph_id)
        else:
            graph_obj = obj.graphs[graph_id]
            plotname = imgdir + "modules_degenr_heatmap_" + \
                suffix + "." + ext
    # Plot DF enrichment based on the graph object:
    statsdf = _plot_df_enrichment_core(
        df=df, col=col, genes=graph_obj.genes,
        module_match=graph_obj.module_match[attr],
        mnames=graph_obj.mnames[attr],
        plotname=plotname, title=title)
    return(statsdf)



def _plot_df_enrichment_core(df, col, genes, module_match, mnames, plotname,
                             title=None):
    """Plot enrichment of modules in dataframe."""
    # Calculate the enrichment:
    if isinstance(df, dict):
        statsdf = None
        pmat, rmat, keys, annot, wr = {}, {}, {}, {}, []
        for dkey in df.keys():
            pmat[dkey], rmat[dkey], keys[dkey], sdf = calc_df_enrichment(
                df[dkey], col, genes, module_match)
            rmat[dkey], annot[dkey] = format_enrichment_matrices(
                rmat[dkey], pmat[dkey])
            wr.append(len(keys[dkey]))
            sdf['dkey'] = dkey
            statsdf = pd.concat([statsdf, sdf])
    else:
        pmat, rmat, keys, statsdf = calc_df_enrichment(
            df, col, genes, module_match)
        rmat, annot = format_enrichment_matrices(rmat, pmat)

    # Plot heatmap of effect size with pvals on top:
    akws = {"ha": "center", "va": "center"}
    if isinstance(df, dict):
        dkey = list(rmat.keys())[0]
        w = 3 + rmat[dkey].shape[1] * len(df)
        h = 1 + 0.2 * (rmat[dkey].shape[0] - 1)
        fig, axs = plt.subplots(1, len(df), figsize=(w, h),
                                gridspec_kw={"width_ratios": wr,
                                             "hspace": 0.025, "wspace": 0.025})
        for i, dkey in enumerate(list(df.keys())):
            yt = False if i > 0 else mnames
            sns.heatmap(
                np.log2(rmat[dkey]), annot=annot[dkey], cmap="RdBu_r",
                yticklabels=yt, xticklabels=keys[dkey],
                ax=axs[i], center=0, cbar=True, fmt="s", annot_kws=akws)
            axs[i].set_aspect(0.5)
            axs[i].set_title(dkey)
    else:
        plt.figure(figsize=(3 + rmat.shape[1],
                            1 + 0.2 * (rmat.shape[0] - 1)))
        ax = plt.gca()
        sns.heatmap(
            np.log2(rmat), annot=annot, cmap="RdBu_r",
            yticklabels=mnames, xticklabels=keys,
            center=0, cbar=True, fmt="s", annot_kws=akws)
        ax = plt.gca()
        ax.set_aspect(0.5)
        if title is not None:
            ax.set_title(title)

    plt.tight_layout()
    plt.savefig(plotname)
    plt.close()
    logging.info(plotname)
    return(statsdf)

