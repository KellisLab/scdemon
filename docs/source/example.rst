.. _example:
Basic example
=============

Basic example using one of the :code:`scanpy` datasets:

.. code-block:: python

    import numpy as np
    import scanpy as sc
    import scdemon as sm
    from scdemon.utils import recipe_full
    from scdemon import plotting as pl

    # Load one of scanpy's datasets and preprocess it
    adata = sc.datasets.pbmc3k()
    recipe_full(adata, preprocess=True, annotate=True)


.. TODO: link to modules objects code documentation
Once we have this dataset, we can create an object that handles any modules and graphs derived from it:

.. code-block:: python

    # Make the modules handling object:
    mod = sm.modules_core(adata, suffix='pbmc_example', svd_k=100)
    mod.setup()


.. TODO: link to graph code documentation
After this, we build a basic gene-gene graph and learn modules on top of it. Any graph we create is an independent :code:`graph` object under the overarching modules object.

.. code-block:: python

    graph_id = 'base'
    mod.make_graph(graph_id)

    # Get the modules and/or print them out:
    mlist = mod.get_modules(graph_id, print_modules=False)
    mod.save_modules(graph_id) # Saves to a tsv file in the current directory


.. TODO: put some of these plots in the documentation
We can also plot these modules in multiple ways:

.. code-block:: python

    # Plot genes on the gene-gene graph and on the gene-level UMAP basis
    pl.plot_genes(mod, graph_id, attr="leiden", show_labels=True, width=16)
    pl.plot_genes(mod, graph_id, basis='umap', attr="leiden", width=16)

    # Plot module expression on the cell-level UMAP basis:
    pl.plot_umap_grid(mod, graph_id)


.. TODO: link gprofiler
We can also use :code:`gprofiler-official` to annotate these modules with functional terms:

.. code-block:: python

    gpres = sm.get_goterms(mod, graph_id)

