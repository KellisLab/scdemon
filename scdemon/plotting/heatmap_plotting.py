#!/usr/bin/python
"""Heatmap plotting: average expression and SVD vs covariates."""
import logging
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
import seaborn as sns


def _plot_svd_corr(cvlist, cv_hratio, cv_mats, cv_ticks, svec, plotname,
                  cbar=False):
    """Plot heatmaps of correlation of svd components and covariates."""
    hfull = np.sum(cv_hratio)
    w = 2.5 * len(svec) / 20 + 1.5 + 0.8 * cbar
    h = 1 * hfull / 6 + 0.1 * len(cvlist)
    fig, axs = plt.subplots(len(cv_mats), 1, figsize=(w, h),
                            gridspec_kw={
                                "height_ratios": cv_hratio,
                                "left": 1.5 / w,
                                "right": 0.99 - 0.8 / w * cbar,
                                "top": 1 - 0.1 / h,
                                "bottom": 0.4 / h})
    # Add heatmaps:
    sfact = np.max(svec) / svec
    for i, covar in enumerate(cvlist):
        cmat = cv_mats[covar].T
        cmat = cmat * sfact[np.newaxis, :]
        sns.heatmap(cmat,
                    yticklabels=cv_ticks[covar],
                    xticklabels=(i == len(cvlist) - 1),
                    cmap="RdBu", ax=axs[i], cbar=cbar, center=0)
        axs[i].set_ylabel(covar, fontsize=12, rotation=0,
                          ha="right", va="center")
    # Save figure:
    plt.tight_layout()
    plt.savefig(plotname)
    plt.close()
    logging.info("Plotted graph to " + plotname)



def plot_svd_corr(obj, cvlist, cbar=False):
    obj.calc_svd_corr(cvlist)
    plotname = obj.imgdir + "svd_corr_" + obj.csuff + ".png"
    obj._calculate_covariate_lengths()
    cv_hratio = [obj.cv_lengths[covar] for covar in cvlist]
    _plot_svd_corr(cvlist=cvlist,
                   cv_hratio=cv_hratio,
                   cv_mats=obj.cv_mats,
                   cv_ticks=obj.cv_ticks,
                   svec=obj.cobj.s,
                   plotname=plotname,
                   cbar=cbar)


# TODO: Keep refactoring suff to graph_id below here:
# TODO: add plotting scale parameters:
# TODO: Clean up - use pre-computed scores!
def plot_heatmap_avgexpr(obj, graph_id, cvlist=None,
                            attr="leiden", cbar=False):
    """Plot heatmap of module average expression in each covariate."""
    plotname = obj.imgdir + "heatmap_" + \
        obj.csuff + "_" + graph_id + "_" + attr + ".png"
    if cvlist is None:
        cvlist = ["celltype", "region", "niareagansc",
                    "cogdx", "Apoe_e4", "msex"]
    obj._calculate_covariate_lengths()
    hratio = [1] + [obj.cv_lengths[covar] for covar in cvlist]
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
    h = 2 * hfull / 6 + 0.1 * len(cvlist)
    fig, axs = plt.subplots(len(cvlist) + 1, 1, figsize=(w, h),
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
    for i, covar in enumerate(cvlist):
        covar_dummy = pd.get_dummies(obj.adata.obs[covar])
        covar_dummy = covar_dummy.to_numpy()
        covar_cols = np.sum(covar_dummy, axis=0)
        tform = covar_dummy / covar_cols[np.newaxis, :]
        avg_lmat = tform.T.dot(scores)
        scaled_lmat = avg_lmat / np.sum(avg_lmat, axis=0)[np.newaxis, :]
        sns.heatmap(scaled_lmat,
                    yticklabels=obj.cv_ticks[covar],
                    xticklabels=(i == len(cvlist) - 1),
                    cmap="Blues", ax=axs[i + 1], cbar=cbar)
        axs[i + 1].set_ylabel(covar, fontsize=14)
    # Save heatmap figure:
    plt.tight_layout()
    plt.savefig(plotname)
    plt.close()
    print("Plotted graph to " + plotname)
