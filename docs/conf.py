# Configuration file for the Sphinx documentation builder.
import os

from sphinx.application import Sphinx

import pairinteraction

# -- Project information -----------------------------------------------------

project = "pairinteraction"
copyright = "2017, Pairinteraction Developers"
author = "Pairinteraction Developers"

version = pairinteraction.__version__  # The short X.Y version, use via |version|
release = version  # The full version, including alpha/beta/rc tags, use via |release|

language = "en"


def write_labels(app: Sphinx, event: str) -> None:
    labels = app.env.domaindata["std"]["labels"]
    with open("labels.txt", "w") as f:
        for label, node in labels.items():
            f.write(f"{label}: {node[1]}\n")


# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.extlinks",
    "sphinx.ext.autosectionlabel",
    "nbsphinx",
    "sphinx.ext.inheritance_diagram",
    "sphinx_autodoc_typehints",
]
autosectionlabel_prefix_document = True
autosectionlabel_maxdepth = 4
templates_path = ["_templates"]
exclude_patterns = ["_build", "_doctrees", "Thumbs.db", ".DS_Store"]  # Ignore these source files and folders
source_suffix = ".rst"
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
autodoc_class_signature = "mixed"  # combine class and __init__ doc
autodoc_typehints = "both"
autodoc_type_aliases = {}  # make type aliases nicer


def setup(app: Sphinx) -> None:
    # trick sphinx autodoc to fully document the classes inside pairinteraction.backend.double
    # instead of just saying 'alias of ...'
    all_pi_types = [
        pairinteraction.backend.float,
        pairinteraction.backend.double,
        pairinteraction.backend.complexfloat,
        pairinteraction.backend.complexdouble,
    ]
    for pi in all_pi_types:
        for obj_name in dir(pi):
            obj = getattr(pi, obj_name)
            name = getattr(obj, "__name__", "")
            if any(name.endswith(suffix) for suffix in ["Float", "Double"]):
                obj.__name__ = obj_name
    # add list of all labels for debugging purposes
    app.connect("build-finished", write_labels)


# -- Options for extlinks -------------------------------------------------
repo_slug = os.getenv("GITHUB_REPOSITORY", "pairinteraction/pairinteraction")
extlinks = {
    "github": (f"https://github.com/{repo_slug}/%s", None),
}


# -- Epilog at the end of every source file -------------------------------------------------
if repo_slug == "pairinteraction/pairinteraction":
    pypi_image = "https://img.shields.io/pypi/v/pairinteraction.svg?color=orange"
    pypi_target = "https://pypi.org/project/pairinteraction/"
else:
    pypi_image = "https://img.shields.io/badge/pypi-TestPyPI-orange.svg?style=flat"
    pypi_target = "https://test.pypi.org/project/pairinteraction/"
user, repo = repo_slug.split("/")

rst_epilog = f"""
.. |gh-workflow| image:: https://github.com/{repo_slug}/actions/workflows/python-wheel.yml/badge.svg
           :target: https://github.com/{repo_slug}/actions/workflows/python-wheel.yml
           :alt: Python Wheel
.. |coverage-cpp-ctest| image:: https://img.shields.io/badge/C%2B%2B_coverage-ctest-blue.svg?style=flat
             :target: https://cuddly-adventure-1w1n2vp.pages.github.io/coverage/cpp-ctest/html/index.html
             :alt: C++ Coverage - ctest
.. |coverage-cpp-pytest| image:: https://img.shields.io/badge/C%2B%2B_coverage-pytest-blue.svg?style=flat
             :target: https://cuddly-adventure-1w1n2vp.pages.github.io/coverage/cpp-pytest/html/index.html
             :alt: C++ Coverage - pytest
.. |coverage-python-pytest| image:: https://img.shields.io/badge/Python_coverage-pytest-blue.svg?style=flat
             :target: https://cuddly-adventure-1w1n2vp.pages.github.io/coverage/python-pytest/html/index.html
             :alt: Python Coverage - pytest
.. |pypi| image:: {pypi_image}
          :target: {pypi_target}
          :alt: PyPI Package
.. |arxiv| image:: /images/arXiv-badge.svg
           :target: https://arxiv.org/abs/1612.08053
           :alt: arXiv:1612.08053
.. |license-lgpl| image:: /images/license-lgpl.svg
             :target: https://www.gnu.org/licenses/lgpl-3.0.html
             :alt: License LGPL v3
.. |license-gpl| image:: /images/license-gpl.svg
             :target: https://www.gnu.org/licenses/gpl-3.0.html
             :alt: License GPL v3
"""
