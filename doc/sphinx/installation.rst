.. _Installation:

Installation
============

This guide provides installation instructions for the `pairinteraction` software, catering to users of all levels.
It offers two primary installation methods:

- **Binary Installers:** Ideal for most users, this method employs `pip` for straightforward installation from the 
  Python Package Index, providing access to the software's Python library and graphical interface.
  OS-specific alternatives are also provided.

- **Building from Source:** Aimed at developers and experienced users, this approach involves compiling the software from source, 
  detailing necessary dependencies, and offering steps for either an automated build with `poetry` / `pip` or a manual build.
  The latter allows for additional customization and the execution of different targets, such as generating documentation.

Binary Installers
-----------------

All binary builds come with an easy-to-use graphical user interface for pair potential calculations, taking into
account electric and magnetic fields. In addition, a Python library is provided which can be used to
write your own code and have more fine-grained control over what pairinteraction does. Quantities as matrix elements,
eigenenergies, or excitation dynamics can be calculated. For usage examples
visit the :ref:`tutorials <Tutorials>` section of the documentation or
download the :github:`Python <tree/master/doc/sphinx/examples_python>`
example programs.

We recommend installing pairinteraction from the `Python Package Index`_ by calling

.. code-block:: bat

        pip install pairinteraction

from the command line. This gives you access to the pairinteraction Python library and the graphical user interface.
The graphical user interface can be started by executing

.. code-block:: bat

    start_pairinteraction_gui

Alternatively, if you are mainly interested in the graphical user interface, you can install a package from :github:`GitHub Releases <releases>` (Windows, OS X) or
from `Flatpak`_ (GNU/Linux). For the installation of the Flatpak package, you have to `install Flatpak`_ first and
then run ``flatpak install org.pairinteraction.Pairinteraction`` from the command line.

.. _Python Package Index: https://pypi.org/project/pairinteraction
.. _Flatpak: https://flathub.org/apps/org.pairinteraction.Pairinteraction
.. _install Flatpak: https://flathub.org/setup


Building from Source
--------------------

Advanced users, especially those who want to contribute to the development of pairinteraction, can build the software from source. The source code is available on
:github:`GitHub <>` and can be cloned from the pairinteraction repository using the following command:

.. code-block:: bat

    git clone --single-branch https://github.com/pairinteraction/pairinteraction.git

The ``--single-branch`` flag is not essential but will speed up the download significantly by omitting all other branches except master.

Requirements
^^^^^^^^^^^^

Before compiling the source code, you have to install the following tools and dependencies:

Build tools
    `CMake`_ for running the build system (at least CMake 3.16 is required), `poetry`_ for managing the Python dependencies

Dependencies for the C++ backend
    If you are using GNU/Linux or OS X, complete lists of dependencies can be found in the Dockerfiles that we use for continuous integration.
    The Dockerfiles are located in the :github:`docker branch <tree/docker/docker>` of the pairinteraction repository. If you are using Windows, you can use `VCPKG`_ with :github:`our configuration file <tree/master/vcpkg.json>` to install the dependencies.

Dependencies for the Python library
    All Python dependencies are listed within the :github:`pyproject.toml <tree/master/pyproject.toml>` file. They are installed automatically when you build the Python library using poetry.

.. _cmake: https://cmake.org
.. _poetry: https://python-poetry.org/docs/#installing-with-the-official-installer
.. _VCPKG: https://vcpkg.io

Automatic Build
^^^^^^^^^^^^^^^

.. note::
    If you do not want to modify the source code and just want to use the most recent version of pairinteraction, you can install pairinteraction directly from the :github:`GitHub <>` repository by running
    ``pip install git+https://github.com/pairinteraction/pairinteraction``. Similarly, you can add the most recent version of pairinteraction to a Python project that is managed by poetry by running ``poetry add git+https://github.com/pairinteraction/pairinteraction``.

After cloning the repository and installing the requirements, you can build and install the software into a local virtual Python environment by running the following command within the pairinteraction repository:

.. code-block:: bat

    poetry install

This will call CMake automatically to build the C++ backend, the Python library, and the graphical user interface. The graphical user interface can be started by executing

.. code-block:: bat

    poetry run start_pairinteraction_gui

To use Python library, you have to run your python code in the virtual environment created by poetry. This can be done by running ``poetry run python your_script.py``.
Alternatively, you can build and install the software system-wide by running ``pip install -e .`` from the root directory of the pairinteraction repository.

Tests of the Python library and graphical user interface can be run by executing

.. code-block:: bat

    poetry run pytest

Manual Build
^^^^^^^^^^^^

.. note::
    Advanced examples for the usage of CMake to build the software for various operating systems can be found in the :github:`workflows <tree/master/.github/workflows>` directory of the pairinteraction repository.

If you want to build, e.g., the documentation of pairinteraction or have more control over the build process, you can run the tasks that have been executed by poetry manually. For this, you have to first install the Python dependencies manually:

.. code-block:: bat

    poetry export -f requirements.txt > requirements.txt
    pip install -r requirements.txt

Then you can build the software using CMake:

.. code-block:: bat

    mkdir build
    cd build
    cmake ..
    cmake --build . --config Release

This creates the C++ backend, the Python library, and the graphical user interface. The graphical user interface can be started by executing

.. code-block:: bat

    ./start_pairinteraction_gui

in the build directory.
To use the Python library, you have to extend the Python package search path to accommodate pairinteraction by adding your build directory to ``PYTHONPATH``.
This can be done e.g. by adding the following lines to the top of a Python script:

.. code-block:: python

    import sys
    sys.path.append("/your/path/to/pairinteraction/build")

Running the different build commands manually has the advantage that you can pass additional options to the build system. For example, you can disable the graphical user interface by running CMake with ``cmake -DWITH_GUI=OFF ..``. A full list of build options is provided in the following:

+---------------------+--------------------------------------+---------+
| Option              | Effect                               | Default |
+=====================+======================================+=========+
| ``WITH_BACKEND``    | Build with C++ backend               | ON      |
+---------------------+--------------------------------------+---------+
| ``WITH_PYTHON``     | Build with SWIG Python interface     | ON      |
+---------------------+--------------------------------------+---------+
| ``WITH_GUI``        | Build with Python GUI                | ON      |
+---------------------+--------------------------------------+---------+
| ``WITH_GSL``        | Use the GNU scientific library for   | ON      |
|                     | Whittaker functions [#]_             |         |
+---------------------+--------------------------------------+---------+
| ``WITH_LAPACKE``    | Use BLAS and LAPACK to speed up      | ON      |
|                     | linear algebra                       |         |
+---------------------+--------------------------------------+---------+
| ``WITH_DOC``        | Generate documentation               | OFF     |
+---------------------+--------------------------------------+---------+
| ``WITH_BENCH``      | Compile the benchmarks               | OFF     |
+---------------------+--------------------------------------+---------+
| ``WITH_DMG``        | Generate a DMG file (Mac OS X only)  | OFF     |
+---------------------+--------------------------------------+---------+
| ``WITH_COVERAGE``   | Generate code coverage report        | OFF     |
+---------------------+--------------------------------------+---------+
| ``WITH_LTO``        | Build with link-time optimization    | OFF     |
+---------------------+--------------------------------------+---------+
| ``WITH_CLANG_TIDY`` | Run Clang-Tidy during compilation    | OFF     |
+---------------------+--------------------------------------+---------+
| ``WITH_JULIA``      | Build a Julia compatible .so         | OFF     |
+---------------------+--------------------------------------+---------+

.. [#] This mode activates the extension for calculating radial wave
       functions using Whittaker functions. If pairinteraction
       is built in this mode, any derived work has to be licensed under
       GPL v3, because of the GSL being distributed under GPL.

Moreover, executing the commands manually allows for running additional targets.
For example, you can use the ``doc`` target to build a documentation by executing ``cmake --build . --target doc``.
In contrast, if you use poetry to build the software, only the default target for building the library is executed.
In the following, a list of all available targets is provided.
Note that some targets require specific build options to be enabled in addition to the default options.

+--------------+-------------------------------------------+----------------------+
| Target       | Task                                      | Requirement          |
+==============+===========================================+======================+
| ``all``      | Build the software (default target)       |                      |
+--------------+-------------------------------------------+----------------------+
| ``test``     | Run the test suite, including C++ tests   |                      |
|              | that are not run by pytest                |                      |
+--------------+-------------------------------------------+----------------------+
| ``bench``    | Run the benchmark suite                   | ``WITH_BENCH=ON``    |
+--------------+-------------------------------------------+----------------------+
| ``doxygen``  | Build the Doxygen documentation           | ``WITH_DOC=ON``      |
|              | in ``doc/doxygen``                        |                      |
+--------------+-------------------------------------------+----------------------+
| ``sphinx``   | Build the Sphinx documentation            | ``WITH_DOC=ON``      |
|              | in ``doc/sphinx``                         |                      |
+--------------+-------------------------------------------+----------------------+
| ``doc``      | Synonym to make both, ``doxygen`` and     | ``WITH_DOC=ON``      |
|              | ``sphinx`` documentation                  |                      |
+--------------+-------------------------------------------+----------------------+
| ``livehtml`` | Build the Sphinx documentation and view   | ``WITH_DOC=ON``      |
|              | it in a web server                        |                      |
+--------------+-------------------------------------------+----------------------+
| ``win32``    | Create a package for Windows              |                      |
+--------------+-------------------------------------------+----------------------+
| ``package``  | Create a packages for GNU/Linux           |                      |
+--------------+-------------------------------------------+----------------------+
| ``package``  | Create a packages for OS X                | ``WITH_DMG=ON``      |
+--------------+-------------------------------------------+----------------------+
| ``license``  | Add the license to a package for OS X     | ``WITH_DMG=ON``      |
+--------------+-------------------------------------------+----------------------+

In addition, a number of options are typically available for the native build tool that is called by CMake.
For example, on GNU/Linux and OS X, you can pass the ``-j num_jobs`` option to the native build tool to enable parallel compilation,
where ``num_jobs`` specifies the maximal number of jobs that will be run. Setting ``num_jobs`` to the number of available
processors can speed up the compilation process significantly.

.. code-block:: bash

    cmake --build . --config Release -- -j 8
