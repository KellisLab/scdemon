#!/usr/bin/env python3
"""Utility scripts for pruning graphs and adjacency matrices."""
import logging
import numpy as np
from scipy import sparse


def prune_graph_degree(graph, degree=0):
    """Prune igraph graph nodes by degree."""
    logging.info("Removing igraph nodes with degree <= " + str(degree))
    import igraph
    todel = [i for i, x in enumerate(graph.degree()) if x <= degree]
    graph.delete_vertices(todel)
    return(graph)


def prune_graph_components(graph, min_size=4):
    """Prune igraph graph components by size."""
    logging.info("Removing graph components smaller than " + str(min_size))
    import igraph
    memb = graph.clusters().membership  # Connected components
    u, c = np.unique(memb, return_counts=True)  # Sizes of components
    comp = set(u[c < min_size])  # Components to remove
    to_delete = np.array([i for i, x in enumerate(memb) if x in comp])
    graph.delete_vertices(to_delete)  # Remove nodes in components
    return(graph)


def prune_adj_degree(aw, degree_cutoff):
    """Remove nodes with degree of `degree_cutoff` or lower."""
    logging.debug("Removing nodes with degree <=" + str(degree_cutoff))
    aw = aw.tocsr()
    rind = np.array((np.sum(aw > 0, axis=1) > degree_cutoff).T)[0]
    cind = np.array((np.sum(aw > 0, axis=0) > degree_cutoff))[0]
    ind = rind + cind
    aw = aw[ind, :]
    aw = aw[:, ind]
    return (aw, ind)


# External / auxiliary functions for graph operations:
def prune_adj_scale(aw, scale=0.90):
    """Prune adjacency matrix by a relative scaling value."""
    logging.debug("Removing edges below " + str(round(scale*100, 2)) +
                 "\\% of max edge for each node.")
    # TODO: handle 0s for max:
    awm1 = np.max(aw, axis=1).T.toarray()[0]
    awm2 = np.max(aw, axis=0).toarray()[0]
    awm = np.vstack((awm1, awm2))
    awm = np.max(awm, axis=0)
    aw = aw.tocoo()
    scl = np.vstack((awm[aw.row], awm[aw.col]))
    scl = np.max(scl, axis=0)
    pct_data = aw.data / scl
    kind = pct_data > scale
    aw = sparse.coo_matrix(
        (aw.data[kind], (aw.row[kind], aw.col[kind])), shape=aw.shape)
    return aw


def prune_adj_knn(aw, k=50, twodir=True, row_only=False):
    """Prune adjacency to maximum k for each node."""
    logging.debug("Pruning graph to maximum k-edges for each node")
    aw = aw.tocoo()
    ind = np.argsort(-aw.data)  # Sort correlations
    mk = np.zeros(aw.shape[0], int)  # Kept margin
    keepind = np.zeros(len(ind), int)
    # Add indices in order of strength:
    for i in range(len(ind)):
        r = aw.row[i]
        c = aw.col[i]
        if twodir:
            cond = mk[r] < k or mk[c] < k
        else:
            cond = mk[r] < k and mk[c] < k
        if cond:
            mk[r] += 1
            if twodir or not row_only:
                mk[c] += 1
            keepind[i] = 1
    # Remove edges:
    kind = ind[keepind == 1]
    aw = sparse.coo_matrix((aw.data[kind], (aw.row[kind], aw.col[kind])),
                           shape=aw.shape)
    return aw

