Managing the database
=====================

The tool allows for

* changing the format of the database downloaded to your computer
* replacing the local database in "data/database" by a shrank version of the downloaded database

The tool is best installed using the uv package manager and meant to be executed from the root directory of the PairInteraction repository.

Installation
------------

After installing PairInteraction itself, run the following command from the root directory of the PairInteraction repository to install the tool:

```bash
uv sync --group tools --inexact
```

Usage
-----

Change the format of the database downloaded to your computer:
```bash
uv run --no-project dbmanager_optimize --compression {UNCOMPRESSED,SNAPPY,ZSTD}
```

To update the local database in "data/database", ensure that you have the most recent database downloaded to your computer and use `dbmanager_shrink` to install it as the local database.

```bash
uv run pairinteraction database remove
uv run pairinteraction database download Rb Sr87_mqdt Sr88_mqdt Sr88_sqdt Yb171_mqdt Yb173_mqdt Yb174_mqdt
uv run --no-project dbmanager_shrink --out data/database
```

The local database is shrunk as follows (see `dbmanager_shrink --help` for the corresponding options):

* Only states with `nu` between 50 and 65 and `l_ryd` of at most 5 are kept, together with the matrix elements
  between them. Tests must not use states outside of this range. The wigner table is filtered accordingly, i.e. it
  only contains the quantum numbers `f` that occur in the kept states.
* In addition, a few low-lying states with `nu` below 10 are kept so that lifetimes and transition rates can be
  calculated. For all kept states, all matrix elements to other kept states are kept.
* The few low-lying states only capture the decay partially, so that the calculated lifetimes are far off. The
  corresponding tests are run nevertheless so that all code is executed, but they only check the resulting values
  if the tests are executed with `pytest --database-dir "" --download-missing`.
  This is also done once in the continuous integration.
