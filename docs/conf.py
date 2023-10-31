# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys
import os


sys.path.insert(0, os.path.abspath("."))
print(sys.path)

import pkg_resources

installed_packages = pkg_resources.working_set
installed_packages_list = sorted(
    ["%s==%s" % (i.key, i.version) for i in installed_packages]
)
print(installed_packages_list)
project = "dynamic-community-fba"
copyright = "2023, S.C.M.A. Wijnen, F. Moro"
authors = ["S.C.M.A. Wijnen", "F. Moro"]

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.napoleon", "sphinx.ext.autodoc"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
