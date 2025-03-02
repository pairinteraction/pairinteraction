.. include:: readme_link.md
   :parser: myst_parser.sphinx_


Documentation
=============

User Guide
    :doc:`overview` - Discover what the software is and its core features.

    :doc:`installation/installation` - Whether you prefer binary installers or compiling from source, we've got all the information you need.

    :doc:`tutorials/tutorials` - Covering both the graphical user interface and the Python library,
    these tutorials enable you to leverage pairinteraction for your projects.

Contributor Guide
    :doc:`contribute/ways` - Discover the many ways you can help improve pairinteraction, from contributing to the repository to providing quantum defects.

    :doc:`contribute/repo` - Ready to dive into development? Here's how to set up your environment for pairinteraction development,
    ensuring you have all the tools you need.

References
    :doc:`API Reference <modules>` - Documentation of classes and functions of pairinteraction's Python library.

Utility Tools
    :doc:`utility/mqdt` - Learn how to calculate states and matrix elements using multi-channel quantum defect theory with our tool written in Julia.

    :doc:`utility/databases` - Find out how to make states and matrix elements available to pairinteraction.

.. toctree::
    :maxdepth: 2
    :caption: User Guide
    :hidden:

    overview.rst
    installation/installation.rst
    tutorials/tutorials.rst
    modules.rst
    references/references.rst

.. toctree::
    :maxdepth: 2
    :caption: Contributor Guide
    :hidden:

    contribute/ways.rst
    contribute/repo.rst
    utility/databases.rst

.. toctree::
    :maxdepth: 2
    :caption: Utility Tools [External Links]
    :hidden:

    mqdt.rst <http://mqdt.pairinteraction.org>
    ryd-numerov <http://numerov.pairinteraction.org>
