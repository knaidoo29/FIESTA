# Configuration file for the Sphinx documentation builder.
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

import vinfo

# -- Project information

project = 'FiEsta'
copyright = '2021-2025, Krishna Naidoo'
author = 'Krishna Naidoo'

version = vinfo.vstr

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx_simplepdf'
]

source_suffix = ['.rst', '.md']

# Napoleon settings
napoleon_numpy_docstring = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']


# -- Options for HTML output

#html_theme = 'sphinx_rtd_theme'
html_theme = 'furo'

html_static_path = ["_static"]
html_theme_options = {
    "light_logo": "FiEsta_logo_small_white.jpg",
    "dark_logo": "FiEsta_logo_small_black.jpg",
}

# -- Options for EPUB output
epub_show_urls = 'footnote'
