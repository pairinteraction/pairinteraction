# noqa: INP001

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


# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.extlinks",
    "nbsphinx",
    "sphinx.ext.inheritance_diagram",
    "sphinx_autodoc_typehints",
    "myst_parser",
]
templates_path = ["_templates"]
exclude_patterns = ["_build", "_doctrees", "Thumbs.db", ".DS_Store"]  # Ignore these source files and folders
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
    </div>
"""


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
