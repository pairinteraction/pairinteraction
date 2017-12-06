Obtaining pairinteraction
=========================

Binary Installer
----------------

Builds are available for GNU/Linux, Mac OS X, and Windows through
`GitHub Releases`_.  The different packages were built on the following
architectures:

.. _GitHub Releases: https://github.com/pairinteraction/pairinteraction/releases

- ``deb`` package: Ubuntu 16.04 amd64
- ``rpm`` package: OpenSUSE Leap x86_64
- Mac OS X ``dmg``: Mac OS X 10.11
- Windows ``exe``: Compiled with Visual Studio 2015

Pairinteraction Python Library
------------------------------

The pairinteraction program comes with a Python 3 library. It can be used for simulating a broad range of two-atom Rydberg systems taking into
account electric and magnetic fields. Quantities as matrix elements, eigenenergies, or excitation dynamics can be calculated. Examples are shown
in the `tutorials`_ section of the documentation. By installing a build from `GitHub Releases`_, the
library gets installed as well. The library requires ``python3``, ``numpy``, ``scipy``, ``matplotlib``,
and ``pyzmq``. After installing these dependencies, the provided Python `example
scripts`_ should work out-of-the-box.

.. _GitHub Releases: https://github.com/pairinteraction/pairinteraction/releases
.. _example scripts: https://github.com/pairinteraction/pairinteraction/tree/master/doc/sphinx/examples_python
.. _tutorials: https://pairinteraction.github.io/pairinteraction/sphinx/html/tutorials.html

For GNU/Linux, all dependencies are installed automatically and the
Python library can be used right away.

For Windows or OS X, we recommend the installation of the Python 3
distribution `Miniconda`_ or `Anaconda`_. Then, the dependencies of the
pairinteraction Python 3 library can be installed by executing the
following command:

.. _Miniconda: https://conda.io/miniconda.html
.. _Anaconda: https://www.anaconda.com/distribution/

.. code-block:: bat

    conda install numpy scipy matplotlib pyzmq

In order to use the pairinteraction Python 3 library with Windows or
OS X, the path to the library has to be added manually to the Python package search path.
Assuming that the pairinteraction software was installed to its default location, this
can be done by adding the following code block to the top of a Python
script:

.. code-block:: python

    import sys
    if sys.platform == "win32": sys.path.append("C:\Program Files\pairinteraction")
    elif sys.platform == "darwin": sys.path.append("/Applications/pairinteraction.app/Contents/Resources")

Alternatively, in case of Windows, you can copy or link the folder ``C:\Program Files\pairinteraction\libpairinteraction`` somewhere into the Python package search path.
In case of OS X, the folder ``/Applications/pairinteraction.app/Contents/Resources/libpairinteraction`` has to be copied or linked.

Building from Source
--------------------

Requirements
^^^^^^^^^^^^

The following tools libraries, including header files, are required
for compiling the source code.

Git
    Git is a version control system used to track changes.

CMake
    The build system is based on CMake.

C++ Backend and Python Interface
""""""""""""""""""""""""""""""""

C++ Compiler
    C++11 capable C++ compiler (e.g., GCC 4.8.1 or later).

Eigen3
    Eigen is a C++ template library for linear algebra: matrices,
    vectors, numerical solvers, and related algorithms.

Sqlite3
   SQLite is a self-contained, high-reliability, embedded,
   full-featured, public-domain, SQL database engine.

Boost
    A library providing advanced C++ features.

GSL
    The GNU Scientific Library (GSL) is a numerical library.

ZeroMQ
    ZeroMQ is a high-performance asynchronous messaging library.

SWIG
    Simplified Wrapper and Interface Generator to generate a Python
    interface from the C++ code.

Python3
    The user interface is provided in Python scripting language.

NumPy and Scipy
    NumPy and SciPy are libraries for the Python programming language, adding
    support for large, multi-dimensional arrays and matrices.

Graphical User Interface
""""""""""""""""""""""""

PyQt5
    Python bindings of the cross-platform GUI toolkit Qt.

Documentation
"""""""""""""

Doxygen
    Doxygen is used as a documentation generator for C++.

Sphinx
    Sphinx is used as a documentation generator for Python.

Build options
^^^^^^^^^^^^^

The build system uses CMake and has some configuration switches.  These are

+-------------------+--------------------------------------+---------+
| Option            | Effect                               | Default |
+===================+======================================+=========+
| ``WITH_BACKEND``  | Build with C++ backend               | ON      |
+-------------------+--------------------------------------+---------+
| ``WITH_GUI``      | Build with Python GUI                | ON      |
+-------------------+--------------------------------------+---------+
| ``WITH_DATABASE`` | Generate the quantum defect database | ON      |
+-------------------+--------------------------------------+---------+
| ``WITH_DOCS``     | Generate documentation               | ON      |
+-------------------+--------------------------------------+---------+
| ``WITH_DMG``      | Generate a DMG file (Mac OS X only)  | OFF     |
+-------------------+--------------------------------------+---------+
| ``WITH_COVERAGE`` | Generate code coverage report        | OFF     |
+-------------------+--------------------------------------+---------+
| ``WITH_ASAN``     | Enable leak detection for unit tests | OFF     |
+-------------------+--------------------------------------+---------+

These options can be passed directly to ``cmake``, i.e.

.. code-block:: bash

    $ cmake -DWITH_GUI=OFF -DWITH_DATABASE=OFF ..

This way we can only build the C++ backend.

Ubuntu 16.04
^^^^^^^^^^^^

Dependencies
""""""""""""

The build system relies on CMake.  To build the Python GUI we need the
PyQT5 toolkit.  The library Eigen3 is header only and thus cross
platform.  Thus you have to install the following packages

.. code-block:: none

    cmake build-essential git libeigen3-dev pyqt5-dev-tools

For the backend we need the following packages

.. code-block:: none

    libboost-all-dev libgsl-dev libsqlite3-dev sqlite3 libzmq3-dev swig python3 python3-dev python3-numpy python3-scipy python3-zmq

The GUI builds with only ``pyqt5-dev-tools`` but to run it we
additionally need

.. code-block:: none

    python3-pint python3-psutil python3-pyqt5.qtsvg

Build Instructions
""""""""""""""""""

To build for GNU/Linux checkout the latest version of pairinteraction
using `git`

.. code-block:: bash

    $ git clone --recursive https://github.com/pairinteraction/pairinteraction.git

.. note::
   Don't forget the ``--recursive`` switch.  Otherwise all the
   submodules will be missing and you won't be able to build
   pairinteraction successfully.

Then proceed with the usual CMake workflow

.. code-block:: bash

    $ cd pairinteraction
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make -j 8

This will build pairinteraction for real and complex matrices.
Afterwards you can start the program from the build directory

.. code-block:: bash

    $ ./pairinteraction

For the GUI to work you need Python3 with the packages ``numpy``,
``scipy``, ``pint``, ``psutil``, and ``pyqt5``.

openSUSE Leap
^^^^^^^^^^^^^

Dependencies
""""""""""""

The build system relies on CMake.  To build the Python GUI we need the
PyQT5 toolkit.  The library Eigen3 is header only and thus cross
platform.  Thus you have to install the following packages

.. code-block:: none

    git cmake eigen3-devel python3-qt5-devel

For the backend we need the following packages

.. code-block:: none

    patterns-openSUSE-devel_C_C++ gcc6-c++ sqlite3 sqlite3-devel boost_1_61-devel gsl-devel zeromq-devel swig python3 python3-devel python3-numpy python3-numpy-devel python3-scipy python3-pyzmq

The GUI builds with only ``python3-qt5-devel`` but to run it we
additionally need

.. code-block:: none

    python3-psutil python3-pip

The package manager ``pip`` is needed to install the ``pint`` package
which we also need

.. code-block:: bash

    $ pip install pint

Build Instructions
""""""""""""""""""

To build for GNU/Linux checkout the latest version of pairinteraction
using ``git``

.. code-block:: bash

    $ git clone --recursive https://github.com/pairinteraction/pairinteraction.git

.. note::
    Don't forget the ``--recursive`` switch.  Otherwise all the submodules
    will be missing and you won't be able to build pairinteraction
    successfully.

Then proceed with the usual CMake workflow

.. code-block:: bash

    $ cd pairinteraction
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make -j 8

This will build pairinteraction for real and complex matrices.
Afterwards you can start the program from the build directory

.. code-block:: bash

    $ ./pairinteraction

For the GUI to work you need Python3 with the packages ``numpy``,
``scipy``, ``pint``, ``psutil``, and ``pyqt5``.

Mac OS X
^^^^^^^^

Dependencies
""""""""""""

Build Instructions
""""""""""""""""""

Code documentation
^^^^^^^^^^^^^^^^^^

To generate the code documentation Doxygen and Sphinx are needed (in
addition to all other build dependencies).  On Ubuntu or Debian it can
be obtained using

.. code-block:: bash

    sudo apt-get install doxygen graphviz python3-sphinx python3-numpydoc

Then checkout the latest version of pairinteraction using `git`

.. code-block:: bash

    $ git clone --recursive https://github.com/pairinteraction/pairinteraction.git

.. note::
    Don't forget the ``--recursive`` switch.  Otherwise all the submodules
    will be missing and you won't be able to build pairinteraction
    successfully.

Then proceed with the usual CMake workflow

.. code-block:: bash

    $ cd pairinteraction
    $ mkdir build
    $ cd build
    $ cmake ..

Instead of calling ``make`` you now call the documentation target

.. code-block:: bash

    $ make doc

.. note::
   You can also build only the Doxygen or the Sphinx documentation by
   building the eponymous targets ``doxygen`` or ``sphinx`` (instead of
   ``doc``).

This will build the documentation in the subdirectory
``doc/doxygen/html`` and ``doc/sphinx/html`` of the build directory.
Open the file ``index.html`` in your browser to see the result.

``make``: Compiling, testing and installing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The command ``make`` is mainly used to compile the source code, but it
can do a number of other things. The generic syntax of the ``make``
command is:

.. code-block:: bash

    $ make [options] [target]

When no target is given, the target ``all`` is used.  The possible
options can be looked up in the ``make`` manual pages.  The following
targets are available:

``all``
    Compiles the complete source code.

``check``
    Runs the testsuite.

``clean``
    Deletes all files that were created during the compilation.

``package``
    On GNU/Linux and Mac OS X this will produce an installable package
    for your platform.

``doxygen``
    Creates the Doxygen code documentation in the ``doc/doxygen``
    subdirectory.

``sphinx``
    Creates the Sphinx code documentation in the ``doc/sphinx``
    subdirectory.

``doc``
    Synonym to make both, ``doxygen`` and ``sphinx``

A number of options are available when calling ``make`` which can be
found in the ``make`` manual pages.  One option we would like to
present here nevertheless which is ``-j num_jobs``, which can be used
for parallel compilation on computers that have more than one CPU or
core.  Here ``num_jobs`` specifies the maximal number of jobs that
will be run.  Setting ``num_jobs`` to the number of available
processors can speed up the compilation process significantly.
