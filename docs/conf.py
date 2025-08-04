# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

# ruff: noqa: INP001

# Configuration file for the Sphinx documentation builder.
import os
from typing import Any

import pairinteraction

# -- Project information -----------------------------------------------------

project = "pairinteraction"
copyright = "2017, Pairinteraction Developers"  # noqa: A001
author = "Pairinteraction Developers"

version = pairinteraction.__version__  # The short X.Y version, use via |version|
release = version  # The full version, including alpha/beta/rc tags, use via |release|

language = "en"

# -- sphinx-polyversion -----------------------------------------------------
if os.getenv("POLYVERSION_DATA"):
    from sphinx_polyversion.api import load  # type: ignore [import-untyped]

    # This adds html_context = {"revisions": [GitRef('main', ...), GitRef('v6.8.9', ...), ...], "current": ...}
    html_context: dict[str, Any] = load()
else:
    html_context = {"current_version": pairinteraction.__version__}

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.extlinks",
    "nbsphinx",
    "sphinx_autodoc_typehints",
    "myst_parser",
    "sphinx_tabs.tabs",
]
templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "_build_polyversion",
    "_doctrees",
    "Thumbs.db",
    ".DS_Store",
]  # Ignore these source files and folders
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"
pygments_style = "sphinx"  # syntax highlighting
todo_include_todos = False


# -- Options for HTML output -------------------------------------------------
html_css_files = ["gallery.css"]

html_theme = "sphinx_rtd_theme"
html_logo = "_static/logo48x48.png"

html_static_path = ["_static"]

htmlhelp_basename = "pairinteractiondoc"


# -- Options for jupyter notebooks -------------------------------------------------
nbsphinx_prolog = """
{% set docname = env.doc2path(env.docname, base=None).split("/")[-1] %}

.. raw:: html

    <style>
      .nbinput .prompt,
      .nboutput .prompt {
        display: none;
      }
    </style>

    <div class="admonition note">
      This page was generated from the Jupyter notebook
      <a class="reference external" href="{{ docname|e }}">{{ docname|e }}</a>.
      Open in
      <a class="reference external" href="https://colab.research.google.com/github/pairinteraction/pairinteraction/blob/master/docs/tutorials/examples_python/{{ docname|e }}">Google Colab</a>.
    </div>
"""  # noqa: E501

# -- Options forautosummary -------------------------------------------
autosummary_ignore_module_all = False


# -- Options for autodoc -------------------------------------------
add_module_names = False  # don't add module names to members
autodoc_class_signature = "separated"  # "mixed": combine class and __init__ doc, "separated": separate them
autodoc_typehints = "both"
autodoc_default_options: dict[str, Any] = {
    "members": True,
    "member-order": "bysource",
    "inherited-members": True,  # include inherited members
    "undoc-members": True,  # include members without docstrings
    "class-doc-from": "class",
}


# -- Options for extlinks -------------------------------------------------
repo_slug = os.getenv("GITHUB_REPOSITORY", "pairinteraction/pairinteraction")
extlinks = {
    "github": (f"https://github.com/{repo_slug}/%s", None),
}


# -- Epilog at the end of every source file -------------------------------------------------
rst_epilog = """
"""
