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

    # Minimal:
    conda create -n scdemonpy -c conda-forge numpy pandas scipy igraph umap-learn leidenalg scanpy seaborn matplotlib

    # If you want to run GO enrichment:
    conda install -n scdemonpy -c bioconda gprofiler-official


The full list of dependencies is as follows:

.. code-block:: console

    # Core:
    numpy, pandas, scipy, igraph, umap-learn, 
    leidenalg, scanpy, seaborn, matplotlib

    # Optional:
    anndata, gprofiler-official, adjustText
