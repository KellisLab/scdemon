#!/usr/bin/python
"""Utility scripts for community detection."""
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

