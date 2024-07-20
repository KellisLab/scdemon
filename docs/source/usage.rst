Usage
=====

.. _installation:
Installation
------------

Install the python package with pip directly from github as: 

.. code-block:: bash

    pip install git+https://github.com/KellisLab/scdemon


For a fresh install, or if version number hasn't changed, first run: `pip uninstall scdemon`

Alternatively, install by cloning the directory as:

.. code-block:: bash

    git clone git@github.com:KellisLab/scdemon.git
    pip install ./scdemon


A minimal conda environment for this package includes:

.. code-block:: console

   conda create -n scdemonpy conda-forge::anndata conda-forge::tqdm conda-forge::pip conda-forge::igraph conda-forge::umap-learn conda::scikit-build


The full list of dependencies is as follows:

.. code-block:: console

    numpy, pandas,
    scipy, tqdm,
    igraph, umap-learn,
    leidenalg, scikit-learn,
    scanpy, anndata,
    seaborn, matplotlib,
    gprofiler-official,
    igraph, adjustText,
    numba
