#!/usr/bin/env python3
"""Utility scripts for community detection."""
import logging
import numpy as np
import leidenalg as la


def compute_leiden_partition(graph,
                             resolution=None,
                             partition_type=None,
                             use_weights=False,
                             n_iterations=-1,
                             random_state=1):
    # Set up the leidenalg arguments:
    partition_kwargs = {}
    partition_kwargs["n_iterations"] = n_iterations
    partition_kwargs["seed"] = random_state

    if partition_type is None:
        if resolution is None:
            partition_type = la.ModularityVertexPartition
        else:
            partition_type = la.RBConfigurationVertexPartition
            partition_kwargs["resolution_parameter"] = resolution

    if use_weights:
        partition_kwargs["weights"] = np.array(
            graph.es["weight"]).astype(np.float64)

    # Run leidenalg to get the partition:
    partition = la.find_partition(
        graph, partition_type, **partition_kwargs)
    return(partition)


def calculate_gene_modules_leiden(graph, resolution=None,
                                    partition_type=None, use_weights=False,
                                    n_iterations=-1, random_state=1):
    """Calculate modules from gene-gene graph using graph clustering."""
    # Create a partition using leiden:
    partition = compute_leiden_partition(
        graph=graph,
        resolution=resolution,
        partition_type=partition_type,
        use_weights=use_weights,
        n_iterations=n_iterations,
        random_state=random_state)
    assign = get_modules_from_partition(graph, partition)
    return(assign)


def get_modules_from_partition(graph, partition):
    # TODO: fix for overlapping clusters
    # NOTE: works for unique assignment but not for overlapping factors
    # Initialize and assign each node to a cluster:
    assign = (np.zeros(len(graph.vs)) - 1).astype(int)
    part_list = [np.array(x) for x in list(partition)]
    nclust = len(part_list)
    for i in range(nclust):
        assign[part_list[i]] = i
    logging.info("Found " + str(nclust) + " clusters")
    return(assign)


def set_module_colors(assign, palette):
    # Assign colors to each cluster:
    # TODO: Fix so it works for overlapping modules too
    nclust = np.max(assign) + 1
    ncols = len(palette)
    rep_col = int(np.ceil((nclust + 2) / (ncols * 1.0)))
    part_cols = np.array((palette * rep_col)[1:(1 + nclust)])
    colors = part_cols[assign]
    return(colors)
