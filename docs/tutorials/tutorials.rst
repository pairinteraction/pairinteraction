.. _tutorials:

Tutorials
=========

The PairInteraction software comes with a graphical user interface (GUI) and a Python API. The GUI provides a quick way
to calculate Rydberg interaction potentials, lifetimes, and energy shifts in the presence of external fields. The Python
API is designed to automate calculations and to have more fine-grained control.

.. _tutorial-api:

Tutorials - Python API
----------------------

**Introduction**

The following Jupyter notebooks introduce the Python API. The other tutorials build on this introduction.

.. nbgallery::

    examples_python/quick_start
    examples_python/concepts

**Examples Demonstrating Basic Usage**

We show exemplary applications of the Python API, ranging from the calculation of single-atom properties such as
lifetimes to calculations of pair potentials.

.. nbgallery::

    examples_python/stark_map
    examples_python/pair_potentials
    examples_python/lifetimes
    examples_python/c3_c6_coefficients
    examples_python/state_atom_object

**Examples Demonstrating Basic Usage**

Some examples showcasing the usage of custom green tensors.

.. nbgallery::

    examples_python/green_tensor/green_tensor_example_cavity
    examples_python/green_tensor/pair_potential_with_green_tensor
    examples_python/green_tensor/c6_near_surfaces.ipynb

**Examples From Publications**

We show how PairInteraction's Python API can be applied to solve complex problems and reproduce results from literature.

.. nbgallery::

    examples_python/effective_hamiltonian
    examples_python/pair_potential_efield_sensitivity
    examples_python/atom_ion_interaction
    examples_python/mqdt

.. _tutorial-gui:

Tutorials - Graphical User Interface
------------------------------------

The graphical user interface (GUI) of PairInteraction can be started by executing

.. code-block:: bash

    pairinteraction gui

from the command line. Below we show a few examples of how one can use the GUI.

.. raw:: html
    :file: examples_gui/slideshow.html
