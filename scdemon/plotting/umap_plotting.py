#!/usr/bin/env python3
"""UMAP plotting - grid of modules."""
import logging
import numpy as np

from matplotlib import pyplot as plt


def _plot_cell_umap(umat, c, plotname=None, ax=None, title=None,
                   s=1, width=8, axlab=False, cmap='viridis'):
    """Plot single umap with given scores as colors."""
    mx = np.max(umat, axis=0)
    mn = np.min(umat, axis=0)
    mid = (mx + mn) / 2
    if ax is None:
        # Define plotting range:
        rn = mx - mn
        height = width * rn[1] / rn[0]
        plt.figure(figsize=(width, height))
        ax = plt.gca()
    if title is not None:
        import textwrap # NOTE: Remove dependency
        tw = textwrap.fill(title, 18)
        ax.text(mid[0], mid[1], tw, fontdict={"fontsize": 7},
                ha='center', va='center')
    ax.set_facecolor("white")
    ax.scatter(umat[:, 0], umat[:, 1], s=s, c=c, marker=".",
               edgecolors="none", cmap=plt.get_cmap(cmap))
    # Remove tick labels:
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    if axlab:
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
    if plotname is not None:
        plt.tight_layout()
        plt.savefig(plotname, dpi=350, bbox_inches="tight")
        plt.close()


def _plot_umap_grid(umat, scores, titles, plotname=None,
                   ind=None, sel=None, width=2, s=0.5):
    """Plot modules scores on UMAP as a grid."""
    nplot = scores.shape[1] if sel is None else len(sel)
    # Determine grid shape: NR/NC
    nr = int(np.round(np.sqrt(nplot) * 0.8))
    nc = int(np.ceil(nplot / (nr * 1.0)))
    # Select specific cells:
    if ind is not None:
        umat = umat[ind, :]
        scores = scores[ind, :]
    # Define plotting range so each plot preserves aspect:
    mx = np.max(umat, axis=0)
    mn = np.min(umat, axis=0)
    rn = mx - mn
    height = width * rn[1] / rn[0]  # Size of each plot
    # Make subplots:
    fig, axs = plt.subplots(nrows=nr, ncols=nc,
                            gridspec_kw={"hspace": 0.01, "wspace": 0.01},
                            figsize=(width * nc, height * nr))
    for i in range(nr):
        for j in range(nc):
            k = i * nc + j
            if nr == 1:
                ax = axs[j]
            elif nc == 1:
                ax = axs[i]
            else:
                ax = axs[i, j]
            if k < nplot:
                k = k if sel is None else sel[k]
                title = None if titles is None else titles[k]
                _plot_cell_umap(umat=umat, c=scores[:, k],
                                title=title, s=s, ax=ax)
            else:
                ax.set_facecolor("white")
                ax.axes.get_xaxis().set_visible(False)
                ax.axes.get_yaxis().set_visible(False)
    # Save figure:
    plt.tight_layout()
    plt.savefig(plotname, dpi=350, bbox_inches="tight")
    plt.close()
    logging.info("Plotted grid of modules on UMAP to: " + plotname)


def plot_umap_grid(obj, graph_id, attr='leiden', plotname=None,
                   ind=None, sel=None, width=2, s=0.5, imgdir='./'):
    """\
        Plot module average expression on the UMAP in a grid

        Parameters
        ----------
        obj : gene_graph | modules
            Object (``gene_graph`` or ``modules``) with modules to plot
        graph_id : str
            Name of graph to work with
        attr : str
            Name of modules in the graph
        plotname
            Name for output file, overriding default naming scheme
        ind : np.array
            Which subset of cells to plot, (default ``None`` is plot all cells)
        sel : list | np.array
            Which set of modules to plot
        width : float
            Width (and height) of each subplot in the grid
        s : float
            Size of points
        imgdir : str
            Directory for images
    """
    if plotname is None:
        plotname = imgdir + "module_umap_grid_" + \
            obj.suffix + "_" + graph_id + ".png"
    # Plot umap grid:
    _plot_umap_grid(umat=obj.adata.obsm["X_umap"],
                    scores=obj.graphs[graph_id].scores[attr],
                    titles=obj.graphs[graph_id].mnames[attr],
                    plotname=plotname,
                    ind=ind, sel=sel, width=width, s=s)

