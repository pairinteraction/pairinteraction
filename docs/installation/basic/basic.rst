.. _basic:

Basic Installation
==================

Installing from Python Package Index
------------------------------------

For an automatic installation, you first need to install `pip`_ in a virtual :ref:`python environment<python_setup>`.

We recommend installing pairinteraction from the `Python Package Index`_ by calling

.. code-block:: bash

    pip install pairinteraction

from the command line. This gives you access to the pairinteraction Python library and the graphical user interface.
The graphical user interface can be started by executing

.. code-block:: bash

    start_pairinteraction_gui

Binary Installers
-----------------

Alternatively, if you are only interested in the graphical user interface, you can download an installer for pairinteraction from :github:`GitHub Releases <releases>` (Windows, OS X) or use the `Flatpak`_ package manager (GNU/Linux). For the installation of the Flatpak package, you have to `install Flatpak`_ first and then run ``flatpak install org.pairinteraction.Pairinteraction`` from the command line.

.. _Python Package Index: https://pypi.org/project/pairinteraction
.. _Flatpak: https://flathub.org/apps/org.pairinteraction.Pairinteraction
.. _install Flatpak: https://flathub.org/setup
.. _pip: https://pypi.org/project/pip/
