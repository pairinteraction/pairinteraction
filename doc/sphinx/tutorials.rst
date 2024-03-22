.. _Tutorials:

Tutorials
=========

The pairinteraction software comes with a graphical user interface (GUI) and a Python API.
The GUI provides a quick way to calculate Rybderg interaction potentials and energy shifts in the presence of external fields.
The Python API is designed to automatize calculations and to have more fine-grained control.

Graphical User Interface
------------------------

.. rubric:: Introduction

Here we introduce how the GUI of pairinteraction can be used to calculate interactions between Rydberg atoms and energy shifts in the presence of applied fields.

.. toctree::
   :maxdepth: 1

   examples_gui/introduction/gui_introduction

.. rubric:: Applications

The following tutorial provides an more sophisticated example on how the GUI can been used for research.

.. toctree::
   :maxdepth: 1

   examples_gui/macrodimers/macrodimers

Python API
----------

.. rubric:: Introduction

Here we show the usage of the API. The first tutorial serves as a quick start guide. The second tutorial provides an in-depth introduction, covering
the basic concepts and topics such as the calculation of C3/C6 coefficients, non-perturbative calculations, and effective Hamiltonians.

.. nbgallery::

   examples_python/quick_start
   examples_python/introduction


.. rubric:: Applications

The following jupyter notebooks show how the pairinteraction Python API can be applied to solve complex problems and reproduce results from literature.

.. nbgallery::

   examples_python/matrix_elements
   examples_python/wavefunctions
   examples_python/comparison_to_saffman_fig13
   examples_python/pair_potential_efield_sensitivity
   examples_python/vdw_near_surface
   examples_python/pair_potential_near_surface
   examples_python/atom_ion_interaction
