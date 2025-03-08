.. include:: readme_link.md
   :parser: myst_parser.sphinx_


Documentation
=============

User Guide
    :doc:`overview` - Discover what the software is and its core features.

    :doc:`installation/installation` - Whether you prefer binary installers or compiling from source, we've got all the information you need.

    :doc:`tutorials/tutorials` - Covering both the graphical user interface and the Python library,
    these tutorials enable you to leverage pairinteraction for your projects.

    :doc:`modules` - Documentation of classes and functions of pairinteraction's Python library.

    :doc:`references/references` - A list of references you can cite when using pairinteraction in your research.

Contributor Guide
    :doc:`contribute/developer` - Information about the pairinteraction library for developers and advanced users who want to contribute to pairinteraction.

    :doc:`contribute/ways` - Discover the many ways you can help improve pairinteraction, from contributing to the repository to providing quantum defects.

    :doc:`contribute/repo` - Ready to dive into development? Here's how to set up your environment for pairinteraction development,
    ensuring you have all the tools you need.

    :doc:`contribute/databases` - Find out how to make states and matrix elements available to pairinteraction.

Utility Tools [External Links]
    `mqdt.jl`_ - Learn how to calculate states and matrix elements using multi-channel quantum defect theory with our tool written in Julia.

    `ryd-numerov`_ - Learn how to calculate states and matrix elements using single quantum defect theory and with our tool written in Python.


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

    contribute/developer.rst
    contribute/ways.rst
    contribute/repo.rst
    contribute/databases.rst

.. toctree::
    :maxdepth: 2
    :caption: Utility Tools [External Links]
    :hidden:

    MQDT.jl <http://mqdt.pairinteraction.org>
    Rydberg Numerov <http://numerov.pairinteraction.org>


.. _ryd-numerov: http://numerov.pairinteraction.org
.. _mqdt.jl: http://mqdt.pairinteraction.org
