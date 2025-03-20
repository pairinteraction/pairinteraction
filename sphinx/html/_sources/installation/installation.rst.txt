.. _installation:

Installation
============


Installing from the Python Package Index (PyPI)
-----------------------------------------------

For users who simply want to use the python interface or the graphical user interface,
we recommend installing pairinteraction from the `Python Package Index`_ by calling

.. code-block:: bash

    pip install --only-binary pairinteraction --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ pairinteraction

from the command line if you are using the `pip`_ package manager.
Alternatively, you can run ``uv pip install [...]``,
if you are using the `uv`_ package manager.

This gives you access to the pairinteraction Python library.
Before starting to use the pairinteraction library, we recommend running the following commands in the terminal
to ensure that the package is correctly installed:

.. code-block:: bash

    # print the current version of the pairinteraction library
    pairinteraction --version

    # run a module test to check if everything is working correctly
    pairinteraction test

    # download the databases for the species you want to use (e.g. Rb and Yb174_mqdt)
    pairinteraction download Rb Yb174_mqdt

For examples on how to use the pairinteraction package, visit the :ref:`tutorial-api` section of the documentation.

In addition, a graphical user interface of pairinteraction is installed and can be started by executing

.. code-block:: bash

    pairinteraction gui

from the command line.
This allows you to do some common and basic calculations without writing any code.
For more information on how to use the graphical user interface, visit the :ref:`tutorial-gui` section of the documentation.



.. _Python Package Index: https://pypi.org/project/pairinteraction
.. _pip: https://pypi.org/project/pip/
.. _uv: https://docs.astral.sh/uv/

.. Binary Installers
.. -----------------

..
    Alternatively, if you are only interested in the graphical user interface,
    you can download an installer for pairinteraction from :github:`GitHub Releases <releases>` (Windows, OS X)
    or use the `Flatpak`_ package manager (GNU/Linux).
    For the installation of the Flatpak package, you have to `install Flatpak`_ first and then
    run ``flatpak install org.pairinteraction.Pairinteraction`` from the command line.

..
    .. _Flatpak: https://flathub.org/apps/org.pairinteraction.Pairinteraction

..
    .. _install Flatpak: https://flathub.org/setup

Advanced Installation
---------------------

For developers and experienced users who want to adjust the source code to their own needs, or want :ref:`to contribute to the package <repository>`
we recommend building the software from source.
Some instructions on how to do this can be found in the :ref:`advanced`.


.. toctree::
    :hidden:

    advanced/advanced.rst
