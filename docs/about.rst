About
=====

What are Rydberg atoms?
-----------------------

Rydberg atoms are atoms with one electron excited to a very high energy level. These atoms have exaggerated properties,
such as a large size, long lifetimes, and strong interactions with other Rydberg atoms, ions, and electric fields. The
properties of the Rydberg atoms can be well-controlled, and the atoms can be trapped individually in so-called optical
tweezers that can be arranged in arbitrary structures. This makes Rydberg atoms a versatile tool for studying quantum
phenomena in fundamental research and for applications in quantum technologies. For example, their strong sensitivity
towards fields are used for quantum sensors. The strong interaction between Rydberg atoms are applied for quantum
simulators and to build quantum gates for quantum computing. While quantum computing platforms based on Rydberg atoms
are still a rather young technology, they already achieve fidelities on par with other leading quantum computing
platforms, and offer the potential for scalable architectures.

For what is the PairInteraction software?
-----------------------------------------

To use Rydberg atoms, it is necessary to understand the interactions between them and how they can be tuned to fulfill
the needs of the desired application. The PairInteraction software package provides a Python library and graphical user
interface to calculate the interaction potentials between Rydberg atoms as well as single-atom properties such as
lifetimes, energies, and energy shifts by electric and magnetic fields. To achieve accurate results, state-of-the-art
approaches are used (such as multi-channel quantum defect theory for atoms with two valence electrons and Green tensor
approaches for calculating interactions). The software is used by researchers around the world to answer questions such
as:

- How do interactions between Rydberg atoms or with electromagnetic fields depend on ...

  - the quantum numbers of the Rydberg atoms?
  - the direction of applied electric and magnetic fields?
  - the geometry of the Rydberg atom arrangement?

- Can interactions be tuned to be robust against small variations of experimental parameters?
- How do Rydberg atoms interact with ions?
- Is the lifetime of the chosen Rydberg state large enough for the desired application?
- Is there a simple, effective Hamiltonian that can describe my system of Rydberg atoms?

The backend of the software is written in C++ and makes use of high-performance libraries such as Intel's MKL library
and the DuckDB database engine. This enables calculations which have been previously rather tedious/impossible like the
accurate calculation of interactions between divalent atoms with complex level structures or pair potentials in presence
of high electric and magnetic fields.

To get started with the software, visit the user guide of our documentation. It contains :ref:`tutorials <tutorials>`
such as a quick start guide, jupyter notebooks with example applications, and an :ref:`API reference <api_reference>`.

Can I contribute to the software?
---------------------------------

The software is open source and licensed under LGPL v3. We are always happy to receive contributions such as pull
requests, bug reports, and new examples for our documentation. Our contributor guide shows :ref:`how to get started as a
contributor <getting_started_as_a_contributor>`. If you are concerned whether the code quality of your contribution
would be sufficient, do not be afraid. We are happy to help. If you have any questions, please do not hesitate to
contact us.
