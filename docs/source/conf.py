# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'cuRDF'
copyright = '2025, Joseph Hart'
author = 'Joseph Hart'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
]
autodoc_mock_imports = [
    "torch",
    "nvalchemiops",
    "MDAnalysis",
    "ase",
    "matplotlib",
    "numpy",
]

templates_path = ['_templates']
exclude_patterns = []

language = 'Python'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
}
