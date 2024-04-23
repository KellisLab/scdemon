#!/usr/bin/python
"""Functions for merging multiple graphs."""
import numpy as np
import leidenalg as la


def make_graphlist(mod, power_list, min_size=4,
                   keep_all_z=True, filter_covariate=None, **kwargs):
    """Make and process a list of graphs."""
    # TODO: allow this to work with a variety of values (SVD k, power, z, etc)
    # Compute the graphs on the modules object:
    graphs = compute_graphs(mod, power_list, keep_all_z=keep_all_z,
                            filter_covariate=filter_covariate, **kwargs)
    # Merge gene sets to ensure same sets:
    full_nodes = get_intersection_graph_names(mod, graphs)
    # Prune graphs to same gene set:
    graphlist = []
    for x in graphs:
        g = mod.graphs[x].graph.copy()
        prune_graph_to_nodes(g, full_nodes, by_index=False)
        graphlist.append(g)
    # Flag nodes not in large components:
    keep_nodes = flag_graphlist_nodes(graphlist, min_size=min_size)
    # Delete vertices not in keep_nodes:
    graphlist = clean_graphlist_nodes(graphlist, keep_nodes, by_index=True)
    return(graphlist, graphs)


def compute_graphs(mod, power_list, keep_all_z=True,
                   filter_covariate=None, **kwargs):
    """Compute graphs along a list of eigenvalue-power values."""
    for power in power_list:
        graph_id = f'p{power}'
        if graph_id not in mod.graphs.keys():
            mod.make_graph(graph_id, multigraph=False,
                           full_graph_only=True, power=power,
                           filter_covariate=filter_covariate,
                           keep_all_z=keep_all_z, **kwargs)
    return(['p' + str(x) for x in power_list])


def get_intersection_graph_names(mod, graphs):
    """Get the names shared by all graphs."""
    full_nodes = None
    for x in graphs:
        nodes = mod.graphs[x].graph.vs['name']
        if full_nodes is None:
            full_nodes = nodes
        else:
            nodes = set(nodes)
            full_nodes = [x for x in full_nodes if x in nodes]
    return(full_nodes)


def flag_graphlist_nodes(graphlist, min_size=4):
    """Return the list of nodes in a min_size+ sized component in any graph."""
    keep_nodes = []
    for g in graphlist:
        memb = g.clusters().membership  # Connected components
        u, c = np.unique(memb, return_counts=True)  # Sizes of components
        comp = set(u[c < min_size])  # Components to remove
        to_keep = [i for i, x in enumerate(memb) if x not in comp]
        keep_nodes = keep_nodes + to_keep
    return(keep_nodes)


def clean_graphlist_nodes(graphlist, keep_nodes, by_index=False):
    """Remove nodes in a graph not in keep_nodes."""
    if type(graphlist) is list:
        for g in graphlist:
            prune_graph_to_nodes(g, keep_nodes, by_index=by_index)
    else:
        # NOTE: Could catch error if not graph, but ok for this purpose
        prune_graph_to_nodes(graphlist, keep_nodes, by_index=by_index)
    return(graphlist)


def prune_graph_to_nodes(graph, keep_nodes, by_index=False):
    # Note could autodetect if by index or name? Issues when names are #s
    if by_index:
        NV = len(graph.vs['name'])
        keep_nodes = np.unique(np.array(keep_nodes))
        to_delete = [i for i in np.arange(NV) if i not in keep_nodes]
    else:
        nodes = graph.vs['name']
        to_delete = [i for i, x in enumerate(nodes) if x not in keep_nodes]
    # Remove these vertices by index:
    graph.delete_vertices(to_delete)


def partition_graphlist(graphlist, resolution=3,
                        n_iterations=-1, random_state=1):
    """Partition nodes on multiple graphs using leidenalg multiplex."""
    partition_type = la.RBConfigurationVertexPartition
    partition_kwargs = {}
    partition_kwargs["n_iterations"] = n_iterations
    partition_kwargs["seed"] = random_state
    partition_kwargs["resolution_parameter"] = resolution
    memb, _ = la.find_partition_multiplex(
        graphlist, partition_type, **partition_kwargs)
    return(memb)
