Overview of PairInteraction's Architecture
==========================================

This page provides an overview of how the PairInteraction software is structured and how it is working internally. The
information is intended for developers and advanced users who want to understand the software's architecture in more
detail.

Code Structure
--------------

The PairInteraction software consists of three parts. Each part is implemented in a separate directory in the
PairInteraction repository:

- :github:`src/cpp <tree/master/src/cpp>` - A high-performance C++ backend for constructing and diagonalizing
  Hamiltonians of systems of Rydberg atoms. Python bindings are created using nanobind_.
- :github:`src/pairinteraction <tree/master/src/pairinteraction>` - A Python library that wraps the Python bindings,
  providing a Python interface for easily accessing the functionality in a pythonic way. Moreover, the Python library
  adds functionality like the perturbative calculation of dispersion coefficients and effective Hamiltonians. The Python
  classes that wrap the Python bindings are defined inside the corresponding submodules ``pairinteraction.ket``,
  ``pairinteraction.basis``, etc. . This classes are defined for the general case of complex-valued data types (and thus
  also use the complex C++ bindings). This classes can be accessed via the ``pairinteraction`` and
  ``pairinteraction.complex`` namespaces. To also make the C++ classes with only real valued data types available, we
  also define e.g. classes like ``BasisAtomReal`` and make them available via the ``pairinteraction.real`` namespace
  directly via ``pairinteraction.real.BasisAtom``. This allows a user to simply switch between real and complex data
  types by changing the import statement from ``import pairinteraction.complex as pi`` to ``import pairinteraction.real
  as pi``, without changing any other code.
- :github:`src/pairinteraction_gui <tree/master/src/pairinteraction_gui>` - A graphical user interface (GUI) that allows
  users to perform common calculations without writing any code. The GUI is built on top of the Python library.

Control Flow
------------

Calculations with the PairInteraction software typically involve three steps that we describe in the following. Hereby,
the result of one step is used as input for the next step. Note that, regardless of whether you use the Python library
or the GUI, the steps are always the same and eventually performed by the C++ backend:

1. Constructing a Basis
~~~~~~~~~~~~~~~~~~~~~~~

A user can specify a basis by restricting the energies and quantum numbers of the states to be included. If a basis for
a single atom should be constructed, the PairInteraction software uses these restrictions internally to search its
:ref:`database <database>` for states that match the criteria. If a basis for two atoms should be constructed, the
PairInteraction software combines eigenstates of two single-atom Hamiltonians so that the resulting pair states match
the criteria.

2. Constructing a Hamiltonian
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Independent of whether a Hamiltonian for one or two atoms should be constructed, the first step is always to obtain
matrix elements of fundamental operators within a single-atom basis (for example, matrix elements for the dipole
operators d0, d+, and d-). The matrix elements of the operators are received from PairInteraction's :ref:`database
<database>` where they are indexed by the corresponding states. Next, the matrix elements of the fundamental operators
are combined to matrix elements of more complex operators (for example, the operator describing the interaction with an
electric field). For the multipole interaction between two atoms, tensor products of the single-atom operators are
required. To accelerate the calculation of the tensor products, we take into account energy restrictions of the pair
basis, ensuring that only matrix elements are constructed that belong to pair states that are energetically allowed.

3. Diagonalizing the Hamiltonian
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Hamiltonian is diagonalized using Intel's MKL library. If multiple Hamiltonians have been constructed, for example,
for different interatomic distances, they can be diagonalized in parallel by passing a list of Hamiltonians to the
diagonalization routine. The parallelization is done in C++, circumventing the global interpreter lock of Python.

.. _nanobind: https://github.com/wjakob/nanobind
