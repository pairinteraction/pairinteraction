Comparing performance of the latest PairInteraction version to other programs
=============================================================================

The tool allows for executing a benchmark

* comparing against other programs
* comparing different floating point data types

The tool is best installed using the uv package manager and meant to be executed from the root directory of the PairInteraction repository.

Installation
------------

After installing PairInteraction itself, run the following command from the root directory of the PairInteraction repository to install the tool:

```bash
uv sync --group tools --inexact
```

Usage
-----

Execute a benchmark comparing against other programs:

```bash
uv run --no-project benchmarking --out results
```

Execute a benchmark comparing different floating point data types:

```bash
uv run --no-project benchmarking --out results --floats
```
