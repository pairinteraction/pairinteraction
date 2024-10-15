# Configuration file for the Sphinx documentation builder.
import os

import pairinteraction

# -- Project information -----------------------------------------------------

project = "pairinteraction"
copyright = "2017, Pairinteraction Developers"
author = "Pairinteraction Developers"

version = pairinteraction.__version__  # The short X.Y version, use via |version|
release = version  # The full version, including alpha/beta/rc tags, use via |release|

language = "en"


# -- General configuration ---------------------------------------------------
# needs_sphinx = '1.0'

extensions = [
    "sphinxcontrib.mermaid",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.extlinks",
    "nbsphinx",
    "sphinxcontrib.autodoc_pydantic",
    "sphinx.ext.graphviz",
    "sphinx.ext.inheritance_diagram",
    "sphinx_autodoc_typehints",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "_doctrees", "Thumbs.db", ".DS_Store"]  # Ignore these source files and folders
mermaid_version = ""  # We want to use the version from cdnjs
source_suffix = ".rst"
master_doc = "index"
pygments_style = "sphinx"  # syntax highlighting
todo_include_todos = False


# -- Options for HTML output -------------------------------------------------
_html_js_attributes = {
    "integrity": "sha512-3j181LWtFFhf1Y8tix6sEqRuN4e9p6V8dH6J6O/bGh5mPk82EA0Y88UZtmlh9awZnhPQqOeB1ogq/NzExIqwKw==",
    "crossorigin": "anonymous",
    "referrerpolicy": "no-referrer",
}
html_js_files = [("https://cdnjs.cloudflare.com/ajax/libs/mermaid/10.7.0/mermaid.min.js", _html_js_attributes)]
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


# -- Options for autodoc and autosummary -------------------------------------------
autoclass_content = "init"
autodoc_default_options = {
    "private-members": False,
    "inherited-members": False,
    "imported-members": False,
}
autosummary_imported_members = False


# -- Options for sphinxcontrib.autodoc_pydantic -------------------------------------------
autodoc_pydantic_model_show_config_summary = False
autodoc_pydantic_model_show_validator_summary = False
autodoc_pydantic_model_show_validator_members = True
autodoc_pydantic_field_list_validators = True
autodoc_pydantic_model_summary_list_order = "bysource"
autodoc_pydantic_model_member_order = "bysource"
autodoc_pydantic_model_erdantic_figure = True


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
.. |linux| image:: https://github.com/{repo_slug}/actions/workflows/linux.yml/badge.svg
           :target: https://github.com/{repo_slug}/actions/workflows/linux.yml
           :alt: Linux
.. |windows| image:: https://github.com/{repo_slug}/actions/workflows/windows.yml/badge.svg
             :target: https://github.com/{repo_slug}/actions/workflows/windows.yml
             :alt: Windows
.. |macos| image:: https://github.com/{repo_slug}/actions/workflows/macos.yml/badge.svg
           :target: https://github.com/{repo_slug}/actions/workflows/macos.yml
           :alt: macOS
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
