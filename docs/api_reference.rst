.. _api_reference:

API Reference
=============

The pairinteraction package provides a Python interface for easily accessing the functionality in a pythonic way. This
Python interface can be accessed via the ``pairinteraction`` module by importing

.. code-block:: python

    import pairinteraction.real as pi

Alternatively, the same functionality can be accessed using complex data types via ``import pairinteraction.complex as
pi``. The two submodules are completely identical in their functionality, but only differ in the data type they use.

All the available classes, methods and functions are documented below:

.. currentmodule:: pairinteraction._wrapped

**Database**

.. autosummary::
    :toctree: _autosummary/

       Database

**Single Atom**

.. autosummary::
    :toctree: _autosummary/

       KetAtom
       BasisAtom
       SystemAtom

**Pair of Atoms**

.. autosummary::
    :toctree: _autosummary/

       KetPair
       BasisPair
       SystemPair
       GreenTensor

**Convenience Functions**

.. autosummary::
    :toctree: _autosummary/

       diagonalize

.. currentmodule:: pairinteraction

**Perturbative Calculations**

.. autosummary::
    :toctree: _autosummary/

       perturbative.C3
       perturbative.C6
       perturbative.EffectiveSystemPair
