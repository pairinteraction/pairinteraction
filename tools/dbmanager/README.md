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
uv run pairinteraction purge
uv run pairinteraction download Rb Sr87_mqdt Sr88_mqdt Sr88_singlet Sr88_triplet Yb171_mqdt Yb173_mqdt Yb174_mqdt
uv run --no-project dbmanager_shrink --out data/database
```
