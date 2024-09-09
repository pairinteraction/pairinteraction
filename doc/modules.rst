API Reference
=============


The heart of the pairinteraction package is a high-performance library for constructing and diagonalizing Hamiltonians of systems of Rydberg atoms. The library is implemented in C++ and can be accessed via the ``pairinteraction`` Python module. Depending on which data type should be used to store the matrix elements of the Hamiltonian, we can choose between four different versions of the library that only differ in the data type used:

.. autosummary::
    :toctree: _autosummary/modules
    :recursive:

    pairinteraction.backend

.. currentmodule:: pairinteraction

Building on top of the high-performance library, the package provides a number of tools for making it easier to work with the library. These tools are implemented in Python and can be accessed via the following submodules:

.. autosummary::
    :toctree: _autosummary/modules
    :recursive:

    pairinteraction.model
    pairinteraction.preprocessing
    pairinteraction.simulation
