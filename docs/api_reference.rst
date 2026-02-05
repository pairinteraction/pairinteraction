.. _api_reference:

API Reference
=============

The PairInteraction package provides a Python interface for easily accessing the functionality in a pythonic way. This
Python interface can be accessed via the ``pairinteraction`` module by importing

.. code-block:: python

    import pairinteraction as pi

Alternatively, the same functionality can be accessed using real data types via ``import pairinteraction.real as pi``
which can accelerate the calculations if no complex numbers are needed (i.e. no fields in y-direction). The classes and
functionality in both cases are identical, only the data type used internally differs.

All the available classes, methods and functions for ``pairinteraction`` (as well as for ``pairinteraction.real``) are
documented below:

.. currentmodule:: pairinteraction

**Database**

.. autosummary::
    :toctree: _autosummary/

       Database

**Single Atom**

.. autosummary::
    :toctree: _autosummary/

       KetAtom
       StateAtom
       BasisAtom
       SystemAtom

**Pair of Atoms**

.. autosummary::
    :toctree: _autosummary/

       KetPair
       StatePair
       BasisPair
       SystemPair

**Convenience Functions**

.. autosummary::
    :toctree: _autosummary/

       diagonalize

**Perturbative Calculations**

.. autosummary::
    :toctree: _autosummary/

       EffectiveSystemPair
       C3
       C6

**Green Tensors**

.. autosummary::
    :toctree: _autosummary/

       green_tensor.GreenTensorFreeSpace
       green_tensor.GreenTensorSurface
       green_tensor.GreenTensorCavity
       green_tensor.GreenTensorInterpolator
