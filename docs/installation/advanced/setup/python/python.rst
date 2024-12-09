.. _python_setup:

Python Setup
============

In order to install the python library, you will need a python environment where you can install all dependencies, especially if you also want to test changes in a :ref:`manual build process <manual>`, or when :ref:`building the documentation <Documentation>` on your own device. Therefore, you first will need to install python on your system. In addition, you will also need a package management system to create virtual environments, such that you can avoid package conflicts with packages you may need for different software projects.
In general, we recommend using `uv`_ due to its speed, but as many users are familiar with `conda`_, we also describe its usage as well.


uv
--

Creating and activating a virtual environment that is able to use all pip functionality with `uv`_ can be done by running the following commands at the root directory of the pairinteraction repository:

.. code-block:: bash

    uv venv --python=3.9 .venv
    source .venv/bin/activate

For Windows users the second command to activate the virtual environment is different: :code:`.venv/Scripts/activate` .

This will create a folder .venv in the pairinteraction directoy which contains all information on the virtual environment.
In order to install all dependencies to build the package, you should run

.. code-block:: bash

    uv pip install -r .build_requirements.txt

.. note::
    If you want to install all dependencies from the `pyproject.toml` file, you should use the following commands, which are different from the usual pip commands:

    .. code-block:: bash

        uv pip compile pyproject.toml --all-extras > requirements.txt
        uv pip install -r requirements.txt


Conda
-----

Creating and activating a virtual environment that is able to use all pip functionality with `conda`_ can be done by running the following commands:

.. code-block:: bash

    conda create --name venv python=3.9
    conda activate venv
    conda install pip-tools

In order to install all dependencies to build the package, you should run from inside the conda environment and the root directory of the pairinteraction repository

.. code-block:: bash

    pip install -r .build_requirements.txt

.. note::
    If you want to install all dependencies from the `pyproject.toml` file, you should use the following commands:

    .. code-block:: bash

        pip-compile pyproject.toml --all-extras --output-file=requirements.txt
        pip install -r requirements.txt


.. _uv: https://pypi.org/project/uv/
.. _pip: https://pypi.org/project/pip/
.. _conda: https://anaconda.org/anaconda/conda
.. _Graphviz: https://www.graphviz.org/
