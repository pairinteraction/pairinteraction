.. _python_setup:

Python Setup
============

In order to install the python library, it is recommended to create a python virtual environment where you can install all dependencies.
Therefore, you will first need a package management system to create the python environment.
This way, you can avoid package conflicts with packages you may need for different software projects.
In general, we recommend using `uv`_ due to its speed, but as many users are familiar with `conda`_, we also describe its usage as well.


uv
--

To install `uv`_ please follow the instructions on the `uv installation page <https://docs.astral.sh/uv/getting-started/installation/>`_.
No further requirements (like installing python) are needed, as `uv`_ takes care of everything for you.

Creating a virtual environment with `uv`_ can be done via

.. code-block:: bash

    uv venv .venv --python 3.10

This will create a folder `.venv/` which contains all information on the virtual environment.
Note, that you can choose a different python version by simply changing the version number.

To activate the virtual environment, you can run

.. code-block:: bash

    # For Linux and MacOS users
    source .venv/bin/activate

    # For Windows users
    .venv\Scripts\activate

In order to install all build dependencies for pairinteraction, run

.. code-block:: bash

    uv pip install -r .build_requirements.txt

inside the pairinteraction repository.

If in addition you want to install all dependencies from the `pyproject.toml` file, you can use the following commands:

.. code-block:: bash

    uv pip compile pyproject.toml --all-extras > requirements.txt
    uv pip install -r requirements.txt

.. _uv: https://docs.astral.sh/uv/


Conda
-----

To install `conda`_ please follow the instructions on the `conda installation page <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_.

Creating and activating a virtual environment that is able to use all pip functionality with `conda`_ can be done by running the following commands:

.. code-block:: bash

    conda create --name venv python=3.10
    conda activate venv
    conda install pip-tools

In order to install all dependencies to build the package, you can run from inside the conda environment and the root directory of the pairinteraction repository

.. code-block:: bash

    pip install -r .build_requirements.txt

If in addition you want to install all dependencies from the `pyproject.toml` file, you should use the following commands:

.. code-block:: bash

    pip-compile pyproject.toml --all-extras --output-file=requirements.txt
    pip install -r requirements.txt

.. _conda: https://anaconda.org/anaconda/conda
