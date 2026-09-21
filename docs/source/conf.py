# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'Introduction to MPI'
copyright = '2026, National Computational Infrastructure'
author = 'NCI Training'

release = '2026'
version = '2026'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.mathjax',
]

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
