.. _automatic:

Automatic Build
===============

**1. Setup**


For the automatic build, you do not need to build the C++ backend manually, but you can work with python only. In order to so, you first have to :ref:`setup your python environment <python_setup>`.
All Python dependencies are listed within the :github:`pyproject.toml <tree/master/pyproject.toml>` file. They are installed automatically when you build the Python library using `pip`_.

.. note::
    If you do not want to modify the source code and just want to use the most recent version of pairinteraction, you can install pairinteraction directly from the :github:`GitHub <>` repository by running ``pip install git+https://github.com/pairinteraction/pairinteraction``.


If not stated otherwise, all commands should be executed from inside the virtual environment and the root directory of the pairinteraction repository.

.. _pip: https://pypi.org/project/pip/

**2. Basic installation**

After cloning the repository and creating the virtual environment, you can build and install the software by running:

.. code-block:: bash

    pip install -e .[gui]

This will call CMake automatically to build the C++ backend, the Python library, and the graphical user interface. In order for pip to find CMake and VCPKG, make sure you have adjusted your environment variables as discussed in our section on :ref:`setting up your development environment<system_setup>`.
By omitting the ``[gui]`` option, you can build the software without installing the additional dependencies needed for the graphical user interface.
The option ``-e`` installs the software in editable mode, which means that changes to the python source code are directly reflected in the installed package.
If you don't want this, you can omit the ``-e`` option to install the current version of the software into the virtual environment.

The graphical user interface can now be started by executing:

.. code-block:: bash

    start_pairinteraction_gui

To use the Python library within your code, you can simply run your python code from inside the virtual environment.

**3. Testing**

First run

.. code-block:: bash

    pip install .[test]

to install the relevant packages for your python environment.
Tests of the Python library and graphical user interface can be run by executing:

.. code-block:: bash

    pytest

If you created your own test and it got skipped, you can get more information on it by running


.. code-block:: bash

    pytest -rsx


**4. Build Documentation**

For building the documentation, we are using `Sphinx`_. As the build process is highly platform dependent, detailed information is provided on our page about the :ref:`documentation <Documentation>`:

.. _Sphinx: https://www.sphinx-doc.org/en/master/index.html

**5. Advanced installation options**

Advanced options for developers when building the package:

.. code-block:: bash

    pip install --no-build-isolation -Cbuild-dir=build_pip -v -e .

| ``--no-build-isolation``: Avoid re-creations of virtual environments for building the package (to use this you first have to install all build dependencies, which are stored inside ``.build_requirements.txt`` so you can install them via ``uv pip install -r .build_requirements.txt``).
| ``-Cbuild-dir=build``: Specify a build directory and reuse it for faster future builds.
| ``-v``: Make the output more verbose.
| ``-e``: Install the package in editable mode (i.e. changes to the python files inside pairinteraction are immediately effective).

To install all dependencies without building the package, confer the :ref:`python setup <python_setup>`.
