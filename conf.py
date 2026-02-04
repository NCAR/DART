# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

# -- Project information -----------------------------------------------------

project = 'DART'
copyright = '2023, University Corporation for Atmospheric Research'
author = 'Data Assimilation Research Section'

# The full version, including alpha/beta/rc tags
release = '11.20.1'
root_doc = 'index'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
     'sphinx_copybutton',
     'sphinx_rtd_theme',
     'sphinx.ext.autodoc',
     'sphinx.ext.mathjax'
]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'models/gitm/testdata1/*', 'models/template/new_model.rst',
]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_show_sphinx = False
html_logo = 'guide/_static/nsf-ncar-dart.png'
html_theme_options = {
    'logo_only': True,
    'includehidden': False
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['guide/_templates']

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['guide/_static']

html_css_files = ['css/custom.css']

# include references
with open(os.path.join(os.path.dirname(__file__), 'guide/references.rst')) as f:
    references_content = f.read()

rst_prolog = references_content
