API Reference
=============


The heart of the pairinteraction package is a high-performance library for constructing and diagonalizing Hamiltonians of systems of Rydberg atoms. The library is implemented in C++ and can be accessed via the following Python module:

.. autosummary::
    :toctree: _autosummary/modules
    :recursive:

    pairinteraction_backend

.. currentmodule:: pairinteraction

Building on top of pairinteraction_backend high-performance library, the package provides a number of tools for making it easier to work with the library. These tools are implemented in Python and can be accessed via the following modules:

.. autosummary::
    :toctree: _autosummary/modules
    :recursive:

    model
    preprocessing
    simulation
