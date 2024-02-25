.. _Development:

Development
===========

TODOs
-----

- Why are the Note/Warning etc. boxes not ending on my local machine?
- Why is the `mermaid` support not working in the `mermaid` test below on my local machine? The text is displayed, but the diagram is not rendered.
- Remove "setup.cfg.in", "setup.py.in", "requirements.txt"
- Make it possible to use poetry, also for managing development dependencies

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
