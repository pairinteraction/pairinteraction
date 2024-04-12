.. _Development:

Development [TODO should be remove before a new release]
========================================================

TODOs
-----

- Why is the `mermaid` support not working in the `mermaid` test below on my local machine? The text is displayed, but the diagram is not rendered.
- Move logic to find dependencies to the top-level `CMakeLists.txt` file
- Do we need individual packages for python distributions or is flatpack enough?
- Remove "setup.cfg.in", "setup.py.in", "requirements.txt"
- Make it possible to use poetry, also for managing development dependencies. Use poetry to setup a python build environment in the CI via `poetry install --no-root`.
- Use pytest to run python tests, see https://stackoverflow.com/questions/52494425/running-a-pytest-test-from-cmake-where-the-test-and-sources-are-in-different-fol
- Adapt workflows so that the correct directories get packed into the releases
- Rename examples folder in sphinx, e.g. examples_cpp to examples/pairinteraction_backend
- Do we need the subfolders cpp and python below pairinteraction_backend?

Literature
----------

- Angular Momentum Coupling and Rabi Frequencies for Simple Atomic Transitions, https://arxiv.org/pdf/0804.4528.pdf

Test of the `mermaid` Support
-----------------------------

.. mermaid::

   sequenceDiagram
      participant Alice
      participant Bob
      Alice->John: Hello John, how are you?
      loop Healthcheck
          John->John: Fight against hypochondria
      end
      Note right of John: Rational thoughts <br/>prevail...
      John-->Alice: Great!
      John->Bob: How about you?
      Bob-->John: Jolly good!
