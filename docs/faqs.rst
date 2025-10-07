FAQs
====

How to use PairInteraction's command line interface?
----------------------------------------------------

PairInteraction's command line interface (CLI) can run the graphical user interface (GUI), tests, and manage the
database. To see documentation of all available commands, run the following command in the terminal:

.. code-block:: bash

    pairinteraction --help

This requires that the CLI is found by the terminal, otherwise use ``python -m pairinteraction --help``.

How to speed up calculations?
-----------------------------

There are several ways to speed up calculations:

- If your Hamiltonian does not need complex numbers (which is typically the case if no fields are pointing along the
  y-direction), use the ``pairinteraction.real`` module instead of ``pairinteraction``.

  .. code-block:: python

      import pairinteraction.real as pi

- Diagonalize many systems in one call. PairInteraction parallelizes the computation over all available CPU cores.

  .. code-block:: python

      systems = [
          pi.SystemPair(basis).set_distance_vector([0, 0, d], unit="micrometer")
          for d in distances
      ]
      pi.diagonalize(systems)

- Pick the ``lapacke_evr`` diagonalization method and restrict the eigenenergies to a tight energy window.

  .. code-block:: python

      pi.diagonalize(
          systems,
          diagonalizer="lapacke_evr",
          energy_range=(min_energy, max_energy),
          energy_range_unit="GHz",
      )

- If a slightly reduced precision is fine with you, use ``float32`` and increase ``rtol`` to allow for a larger relative
  error in the eigenenergies.

  .. code-block:: python

      pi.diagonalize(
          systems,
          float_type="float32",  # float64 is the default
          rtol=1e-5,  # 1e-6 is the default
      )

- Always start with a strongly restricted basis (typically, it is a good idea to restrict the quantum numbers ``n`` and
  ``l`` for the single-atom basis and the ``energy`` for the two-atom basis). Then, loosen the restrictions until
  calculations have converged.

How to fix database problems?
-----------------------------

PairInteraction uses a database of states and matrix elements to construct Hamiltonians. Because the database tables are
large, they are not included in the PairInteraction installation but are downloaded once they are needed. The following
guide helps to solve problems with the database.

- If you see errors about missing species/tables (e.g., **“No table … found”**) and you are using the Python API, ensure
  that you have allowed PairInteraction to download tables automatically. If you are using the GUI, missing tables are
  always downloaded.

  .. code-block:: python

      pi.Database.initialize_global_database(download_missing=True)

  Alternatively, you can download tables manually with PairInteraction's CLI using the species identifiers

  .. code-block:: bash

      $ pairinteraction download Rb Yb171_mqdt

  or a URL to a table file from https://github.com/pairinteraction/database-sqdt/releases or
  https://github.com/pairinteraction/database-mqdt/releases

  .. code-block:: bash

      $ pairinteraction download https://github.com/pairinteraction/database-sqdt/releases/download/v1.2/Rb_v1.2.zip

- If you see access errors (e.g., **"Rate limit reached ..."**), GitHub might have temporary blocked the download of
  database tables for unauthenticated users for whom very strict `rate limits`_ apply. The rate limits should typically
  reset within one hour. Alternatively, you can create a `GitHub Personal Access Token`_ for authentication. To make
  PairInteraction use the token, set the environment variable ``GITHUB_TOKEN`` to the value of the token.
- If a problem with the database persists, you can try whether purging the cache via

  .. code-block:: bash

      $ pairinteraction purge

  helps. Note that this deletes all cached files and database tables must be re-downloaded. To manually inspect the
  content of the cache, you can obtain its path via ``pairinteraction paths``.

.. _github personal access token: https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens

.. _rate limits: https://docs.github.com/en/rest/using-the-rest-api/rate-limits-for-the-rest-api

How to obtain debug information?
--------------------------------

If a problem occurs and you want to investigate it further, the following steps can be useful. They produce output that
can help tracking down a bug and :github:`opening an issue on GitHub <issues>`.

- Increase the verbosity of the logging output. If you are using the GUI, run it in the debug mode.

  .. code-block:: bash

      pairinteraction --log-level DEBUG gui

  If you are using the Python API, add the following code to the beginning of your program:

  .. code-block:: python

      import logging

      logging.basicConfig(level=logging.DEBUG)

- Use PairInteraction's CLI to run a self test

  .. code-block:: bash

      pairinteraction --log-level DEBUG test
