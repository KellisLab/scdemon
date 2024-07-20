#!/usr/bin/env python3
"""Heatmap plotting: average expression and SVD vs covariates."""
import logging
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
import seaborn as sns


# TODO: standardize imgdir

# Plot given the calculated matrices and relevant attributes
def _plot_svd_corr(covariate_list, covariate_hratio, covariate_matrices,
                   covariate_ticks, svec, plotname, cbar=False):
    """Plot heatmaps of correlation of svd components and covariates."""
    hfull = np.sum(covariate_hratio)
    w = 2.5 * len(svec) / 20 + 1.5 + 0.8 * cbar
    h = 1 * hfull / 6 + 0.1 * len(covariate_list)
    fig, axs = plt.subplots(len(covariate_matrices), 1, figsize=(w, h),
                            gridspec_kw={
                                "height_ratios": covariate_hratio,
                                "left": 1.5 / w,
                                "right": 0.99 - 0.8 / w * cbar,
                                "top": 1 - 0.1 / h,
                                "bottom": 0.4 / h})
    # Add heatmaps:
    sfact = np.max(svec) / svec
    for i, covar in enumerate(covariate_list):
        cmat = covariate_matrices[covar].T
        cmat = cmat * sfact[np.newaxis, :]
        sns.heatmap(cmat,
                    yticklabels=covariate_ticks[covar],
                    xticklabels=(i == len(covariate_list) - 1),
                    cmap="RdBu", ax=axs[i], cbar=cbar, center=0)
        axs[i].set_ylabel(covar, fontsize=12, rotation=0,
                          ha="right", va="center")
    # Save figure:
    plt.tight_layout()
    plt.savefig(plotname)
    plt.close()
    logging.info("Plotted graph to " + plotname)



# Plot from object
def plot_svd_corr(obj, covariate_list, cbar=False, imgdir='./'):
    obj._calculate_covariate_svd_correlation()
    obj._calculate_covariate_lengths()
    plotname = imgdir + "svd_corr_" + obj.suffix + ".png"
    covariate_hratio = [obj.covariate_lengths[covar]
                        for covar in covariate_list]
    _plot_svd_corr(covariate_list=covariate_list,
                   covariate_hratio=covariate_hratio,
                   covariate_matrices=obj.covariate_matrices,
                   covariate_ticks=obj.covariate_ticks,
                   svec=obj.s,
                   plotname=plotname,
                   cbar=cbar)


# TODO: Keep refactoring suff to graph_id below here:
# TODO: add plotting scale parameters:
# TODO: Clean up - use pre-computed scores!
def plot_heatmap_avgexpr(obj, graph_id, covariate_list=None,
                         attr="leiden", cbar=False, imgdir='./'):
    """Plot heatmap of module average expression in each covariate."""
    plotname = imgdir + "heatmap_" + \
        obj.suffix + "_" + graph_id + "_" + attr + ".png"
    if covariate_list is None:
        covariate_list = ["celltype", "region", "niareagansc",
                    "cogdx", "Apoe_e4", "msex"]
    obj._calculate_covariate_lengths()
    hratio = [1] + [obj.covariate_lengths[covar] for covar in covariate_list]
    hfull = np.sum(hratio)
    # For subgraph colors heatmap:
    partition = obj.graphs[graph_id].assign[attr]
    part_cols = obj.graphs[graph_id].colors[attr]
    u, c = np.unique(partition, return_index=True)
    col_mat = u[np.newaxis, :]
    col_cmap = sns.color_palette(part_cols[c])
    # Make heatmap figure:
    scores = obj.graphs[graph_id].scores[attr]
    w = 7 * scores.shape[1] / 20 + 1.2 + 0.8 * cbar
    h = 2 * hfull / 6 + 0.1 * len(covariate_list)
    fig, axs = plt.subplots(len(covariate_list) + 1, 1, figsize=(w, h),
                            gridspec_kw={"height_ratios": hratio,
                                            "hspace": 0.05,
                                            "left": 1.2 / w,
                                            "right": 1 - 0.8 / w,
                                            "top": 1 - 0.1 / h,
                                            "bottom": 0.4 / h})
    # Plot the module colors:
    sns.heatmap(col_mat, yticklabels=False, xticklabels=False,
                cmap=col_cmap, ax=axs[0], cbar=False)
    # NOTE: only works for categorical covariates:
    for i, covar in enumerate(covariate_list):
        covar_dummy = pd.get_dummies(obj.adata.obs[covar])
        covar_dummy = covar_dummy.to_numpy()
        covar_cols = np.sum(covar_dummy, axis=0)
        tform = covar_dummy / covar_cols[np.newaxis, :]
        avg_lmat = tform.T.dot(scores)
        scaled_lmat = avg_lmat / np.sum(avg_lmat, axis=0)[np.newaxis, :]
        sns.heatmap(scaled_lmat,
                    yticklabels=obj.covariate_ticks[covar],
                    xticklabels=(i == len(covariate_list) - 1),
                    cmap="Blues", ax=axs[i + 1], cbar=cbar)
        axs[i + 1].set_ylabel(covar, fontsize=14)
    # Save heatmap figure:
    plt.tight_layout()
    plt.savefig(plotname)
    plt.close()
    print("Plotted graph to " + plotname)
