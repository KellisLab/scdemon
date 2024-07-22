# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import setuptools

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
needs_sphinx = "4.0"  # For param docs

project = 'scdemon'
copyright = '2024, Benjamin James, Carles Boix'
author = 'Benjamin James, Carles Boix'
repository_url = "https://github.com/KellisLab/scdemon"

version = setuptools.version.metadata.version('scdemon')
release = version


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'myst_nb',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',

    "sphinx.ext.intersphinx",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.extlinks",
    # "sphinxcontrib.bibtex",
]

templates_path = ['_templates']
exclude_patterns = []


myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "html_admonition",
]

# Options cribbed from scanpy docs
autosummary_generate = True
autodoc_member_order = "bysource"
autosummary_imported_members = True # For class
# autodoc_default_flags = ['members']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
todo_include_todos = False

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
# NOTE: Scanpydoc theme is ok too -
# (sphinx-book-theme, with readthedocs-sphinx-search)
# html_theme = "scanpydoc"
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
# html_theme_options = { # Meeds sphinx-book-theme or scanpydoc
#     "repository_url": repository_url,
#     "use_repository_button": True,
# }
html_show_sphinx = False
# html_logo = "_static/img/Scanpy_Logo_BrightFG.svg"
html_title = "scdemon"

