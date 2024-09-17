API Reference
=============


The heart of the pairinteraction package is a high-performance library for constructing and diagonalizing Hamiltonians of systems of Rydberg atoms.
The library is implemented in C++ and can be accessed via the ``pairinteraction.backend`` Python module.

.. currentmodule:: pairinteraction

.. autosummary::
    :toctree: _autosummary/modules
    :recursive:

    backend

Building on top of the high-performance library, the package provides a number of tools for making it easier to work with the library.
These tools are implemented in Python and can be accessed via the following submodules (e.g. ``pairinteraction.simulation``):

.. autosummary::
    :toctree: _autosummary/modules
    :recursive:

    model
    preprocessing
    simulation
