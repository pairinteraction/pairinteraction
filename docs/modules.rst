API Reference
=============


The heart of the pairinteraction package is a high-performance library for constructing and diagonalizing Hamiltonians of systems of Rydberg atoms implemented in C++.

The library provides a Python interface for easily accessing the functionality in a pythonic way.
This Python interface can be accessed via the ``pairinteraction`` module by importing

.. code-block:: python

    import pairinteraction.backend.real as pi

Alternatively, the same functionality can be accessed for using complex data type by replacing ``real`` with ``complex``.
These four submodules are completely identical in their functionality, but only differ in the data type they use in the C++ backend.
All the available classes, methods and functions are documented here for the ``real`` data type, but the same classes, methods and functions are available for the complex data types as well.

.. autosummary::
    :toctree: _autosummary/modules
    :recursive:

    pairinteraction.backend.real


.. note::

    For advanced users and developers:
    The Python classes are defined inside ``pairinteraction.backend._wrapped`` and then aliased to the ``pairinteraction.backend.real`` and ``pairinteraction.backend.complex`` namespaces.
    The bare C++ functionality is bound to Python using `nanobind`_, all the original bound C++ classes and functions can be accessed via the ``pairinteraction.backend._backend`` namespace.


.. _nanobind: https://github.com/wjakob/nanobind
