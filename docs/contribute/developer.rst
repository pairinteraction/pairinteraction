Developer Informations
======================

This page contains additional information about the pairinteraction library, which is intended for developers and advanced users.

The heart of the pairinteraction library is a high-performance library for constructing and diagonalizing Hamiltonians of systems of Rydberg atoms implemented in C++.

The library provides a Python interface for easily accessing the functionality in a pythonic way.
The Python classes are defined inside ``pairinteraction._wrapped`` and then aliased to the ``pairinteraction.real`` and ``pairinteraction.complex`` namespaces.
The two submodules ``real`` and ``complex`` are completely identical in their functionality, but only differ in the data type they use in the C++ backend.

The bare C++ functionality is bound to Python using `nanobind`_, all the original bound C++ classes and functions can be accessed via the ``pairinteraction._backend`` namespace.

.. _nanobind: https://github.com/wjakob/nanobind
