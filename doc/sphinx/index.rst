**************************************************
Pairinteraction - A Rydberg Interaction Calculator
**************************************************

|linux| |windows| |macos| |codecov| |benchs| |pypi| |arxiv| |license|

TODO banner picture !!!!!!!!!!!!!!!!!!!!!

TODO extend contributing to repo !!!!!!!!!!!!!!

TODO collapsable reference list for quantum defects !!!!!!!!!!!!!!

TODO use index.rst also as README? !!!!!!!!!!!!!!

The *pairinteraction* software calculates properties of Rydberg systems.
The software consists of a C++ backend, a Python library, and a graphical user interface for pair potential calculations.
For usage examples visit the :ref:`tutorials <Tutorials>` section of the documentation.
Stay tuned by `signing up`_ for the newsletter so whenever there are updates to the software or new publications about pairinteraction we can contact you.
If you have a question that is related to problems, bugs, or suggests an improvement, consider raising an :github:`issue <issues>` on :github:`GitHub <>`.

Binary builds are available through :github:`GitHub Releases <releases>`. 
For using pairinteraction as a Python 3 library, we recommend the installation via pip by calling

.. code-block:: bash

    pip install pairinteraction

If pairinteraction was installed via pip, the graphical user
interface can be started by executing ``start_pairinteraction_gui`` from the command line.

**Please cite as**

    Sebastian Weber, Christoph Tresp, Henri Menke, Alban Urvoy, Ofer Firstenberg, Hans Peter Büchler, Sebastian Hofferberth,
    *Tutorial: Calculation of Rydberg interaction potentials*,
    `J. Phys. B: At. Mol. Opt. Phys. 50, 133001 (2017) <https://doi.org/10.1088/1361-6455/aa743a>`_, `arXiv:1612.08053 <https://arxiv.org/abs/1612.08053>`_ 

The pairinteraction software relies on quantum defects provided by the community.
Please consider citing relevant publications for your atomic species alongside pairinteraction:
`Rb` [add link to copy a bibtex entry], TODO add all species with the name that we use within the software.

.. _signing up: https://goo.gl/forms/4bmz3qeuLjKfRlWJ3
.. _LGPL v3: https://www.gnu.org/licenses/lgpl-3.0.html
.. _GPL v3: https://www.gnu.org/licenses/gpl-3.0.html

Documentation
=============

User Guide
    Use our comprehensive user guide, designed to help you effectively utilize the software.

    :doc:`overview` - Discover what the software is and its core features.

    :doc:`installation` - Whether you prefer binary installers or compiling from source, we've got all the information you need.

    :doc:`tutorials` - Covering both the graphical user interface and the Python library, 
    these tutorials enable you to leverage pairinteraction for your projects.

Contributer Guide
    Interested in contributing to the Pairinteraction project? We welcome contributions from everyone, regardless of your experience level!

    :doc:`ways` - Not sure where to start? Discover the many ways you can contribute to improving pairinteraction, 
    from contributions to the repository to providing quantum defects.

    :doc:`repo` - Ready to dive into development? Here's how to set up your environment for pairinteraction development, 
    ensuring you have all the tools you need.

:doc:`API Reference <modules>`
    Documentation of classes and functions of pairinteraction's Python library.

:doc:`Config Reference <config>`
    Documentation of the configuration file for describing systems of Rydberg atoms.

Utility Tools
    Pairinteraction comes with a suite of utility tools to address specific tasks.

    :doc:`mqdt` - Learn how to calculate states and matrix elements using multi-channel quantum defect theory with our tool written in Julia.

.. toctree::
    :numbered:
    :maxdepth: 2
    :caption: User Guide
    :hidden:

    overview.rst
    installation.rst
    tutorials.rst
    development.rst

.. toctree::
    :numbered:
    :maxdepth: 2
    :caption: Contributer Guide
    :hidden:

    ways.rst
    repo.rst

.. toctree::
    :maxdepth: 2
    :caption: API Reference
    :hidden:

    modules.rst

.. toctree::
    :maxdepth: 2
    :caption: Config Reference
    :hidden:

    config.rst

.. toctree::
    :maxdepth: 2
    :caption: Utility Tools
    :hidden:

    mqdt.rst

Credits
=======

TODO: Update credits. Maintainers, list of contributors.

The pairinteraction software was originally developed at the `5th Institute of Physics`_ and the `Institute for Theoretical Physics III`_ of the University of Stuttgart, Germany.
Currently it is maintained by developers at the `Institute for Theoretical Physics III`_ of the University of Stuttgart in Germany, the `Department of Physics`_ of the University of Otago in New Zealand,
the `Institute of Physics`_ of the University of Rostock in Germany, and the `Department of Physics, Chemistry and Pharmacy`_ of the University of Southern Denmark in Denmark.

.. _5th Institute of Physics: http://www.pi5.uni-stuttgart.de/
.. _Institute for Theoretical Physics III: http://www.itp3.uni-stuttgart.de/
.. _Department of Physics: http://www.otago.ac.nz/physics/index.html
.. _Department of Physics, Chemistry and Pharmacy: http://www.sdu.dk/en/fkf
.. _Institute of Physics: https://www.physik.uni-rostock.de/

License
=======

The pairinteraction library is licensed under the `LGPL v3`_. The extension for calculating
radial wave functions using Whittaker functions and the graphical user interface are licensed under the `GPL v3`_.
The GPL v3 also applies to the combined work and all provided binary builds.