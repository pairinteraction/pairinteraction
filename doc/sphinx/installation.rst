Installation
============

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

The binary builds come with a Python and C++ library. The library can be used
for simulating a broad range of two-atom Rydberg systems taking into
account electric and magnetic fields. Quantities as matrix elements,
eigenenergies, or excitation dynamics can be calculated. For usage examples
visit the `tutorials`_ section of the documentation or download the `Python`_
and `C++`_ example programs.

.. _tutorials: https://pairinteraction.github.io/pairinteraction/sphinx/html/tutorials.html
.. _python: https://github.com/pairinteraction/pairinteraction/tree/master/doc/sphinx/examples_python
.. _C++: https://github.com/pairinteraction/pairinteraction/tree/master/doc/sphinx/examples_cpp

Python Library
^^^^^^^^^^^^^^

For GNU/Linux, all dependencies of the pairinteraction Python 3 library are installed
automatically and the Python library can be used right away.

For Windows or OS X, we recommend the installation of the Python 3
distribution `Miniconda`_ or `Anaconda`_. Then, the dependencies of the
pairinteraction Python 3 library can be installed by executing the
following command:

.. _Miniconda: https://conda.io/miniconda.html
.. _Anaconda: https://www.anaconda.com/distribution/

.. code-block:: bat

    conda install numpy scipy matplotlib

In order to use the pairinteraction Python 3 library with Windows or
OS X, the path containing the library has to be added manually to the Python package search path.
Assuming that the pairinteraction software was installed to its default location, this
can be done by adding the following code block to the top of a Python
script:

.. code-block:: python

    import sys
    if sys.platform == "win32": sys.path.append("C:\Program Files\pairinteraction")
    elif sys.platform == "darwin": sys.path.append("/Applications/pairinteraction.app/Contents/Resources")

C++ Library
^^^^^^^^^^^

For GNU/Linux, all dependencies of the pairinteraction C++ library are installed
automatically and the C++ library can be used right away.

For Windows or OS X, the following C++ libraries and their development headers have to
be installed manually: Sqlite3, Boost, GSL

Building from Source
--------------------

.. warning::
   In the tutorials and examples you will often find instructions to
   extend the ``PYTHONPATH`` to accomodate pairinteraction.  These
   instructions only hold true if you installed pairinteraction from
   the official binary installers.  If you built from source you
   instead have to point the scripts to your build directory or the
   custom location where you installed the program.

Requirements
^^^^^^^^^^^^

The following tools and libraries, including header files, are required
for compiling the source code.

Git
    Git is a version control system used to track changes.

CMake
    The build system is based on CMake.

C++ Backend and Python Interface
""""""""""""""""""""""""""""""""

C++ Compiler
    C++11 capable C++ compiler (e.g., GCC 4.8.1 or later).

Sqlite3
   SQLite is a self-contained, high-reliability, embedded,
   full-featured, public-domain, SQL database engine.

Boost
    A library providing advanced C++ features.

GSL
    The GNU Scientific Library (GSL) is a numerical library.

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

+---------------------+--------------------------------------+---------+
| Option              | Effect                               | Default |
+=====================+======================================+=========+
| ``WITH_BACKEND``    | Build with C++ backend               | ON      |
+---------------------+--------------------------------------+---------+
| ``WITH_PYTHON``     | Build with SWIG Python interface     | ON      |
+---------------------+--------------------------------------+---------+
| ``WITH_GUI``        | Build with Python GUI                | ON      |
+---------------------+--------------------------------------+---------+
| ``WITH_DATABASE``   | Generate the quantum defect database | ON      |
+---------------------+--------------------------------------+---------+
| ``WITH_DOCS``       | Generate documentation               | ON      |
+---------------------+--------------------------------------+---------+
| ``WITH_DMG``        | Generate a DMG file (Mac OS X only)  | OFF     |
+---------------------+--------------------------------------+---------+
| ``WITH_COVERAGE``   | Generate code coverage report        | OFF     |
+---------------------+--------------------------------------+---------+
| ``WITH_CLANG_TIDY`` | Run Clang-Tidy during compilation    | OFF     |
+---------------------+--------------------------------------+---------+

These options can be passed directly to ``cmake``, i.e.

.. code-block:: bash

    $ cmake -DWITH_GUI=OFF -DWITH_DATABASE=OFF ..

This way we can only build the C++ backend.

Ubuntu 16.04
^^^^^^^^^^^^

Dependencies
""""""""""""

The build system relies on CMake.  To build the Python GUI we need the
PyQT5 toolkit.  Thus you have to install the following packages

.. code-block:: none

    cmake build-essential git pyqt5-dev-tools

For the backend we need the following packages

.. code-block:: none

    libboost-all-dev libgsl-dev libsqlite3-dev sqlite3 swig python3 python3-dev python3-numpy python3-scipy

The GUI builds with only ``pyqt5-dev-tools`` but to run it we
additionally need

.. code-block:: none

    python3-pint python3-psutil python3-pyqt5.qtsvg

Build Instructions
""""""""""""""""""

To build for GNU/Linux checkout the latest version of pairinteraction
using `git`

.. code-block:: bash

    $ git clone --single-branch --recursive https://github.com/pairinteraction/pairinteraction.git

.. note::
   Don't forget the ``--recursive`` switch.  Otherwise all the
   submodules will be missing and you won't be able to build
   pairinteraction successfully.  The ``--single-branch`` flag is not
   essential but will speed up the download significantly by omitting
   all other branch except ``master``.

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

For the GUI to work, you need Python3 with the packages ``numpy``,
``scipy``, ``pint``, ``psutil``, and ``pyqt5``.

In order to use the pairinteraction Python 3 library,
you have to add the build directory to the Python package search path. The pairinteraction C++
library can be used right away.

openSUSE Leap
^^^^^^^^^^^^^

Dependencies
""""""""""""

The build system relies on CMake.  To build the Python GUI we need the
PyQT5 toolkit.  Thus you have to install the following packages

.. code-block:: none

    git cmake python3-qt5-devel

For the backend we need the following packages

.. code-block:: none

    patterns-openSUSE-devel_C_C++ gcc6-c++ sqlite3 sqlite3-devel boost_1_61-devel gsl-devel swig python3 python3-devel python3-numpy python3-numpy-devel python3-scipy

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

    $ git clone --single-branch --recursive https://github.com/pairinteraction/pairinteraction.git

.. note::
    Don't forget the ``--recursive`` switch.  Otherwise all the
    submodules will be missing and you won't be able to build
    pairinteraction successfully.  The ``--single-branch`` flag is not
    essential but will speed up the download significantly by omitting
    all other branch except ``master``.

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

For the GUI to work, you need Python3 with the packages ``numpy``,
``scipy``, ``pint``, ``psutil``, and ``pyqt5``.

In order to use the pairinteraction Python 3 library,
you have to add the build directory to the Python package search path. The pairinteraction C++
library can be used right away.

Mac OS X
^^^^^^^^

Dependencies
""""""""""""

The build system relies on CMake. For building the pairinteraction C++ backend,
you have to install (e.g. via homebrew) the following packages

.. code-block:: none

    cmake git gsl swig llvm@3.9

.. note::
    The package llvm contains the Clang C++ compiler. We use this compiler as it
    supports OpenMP with OS X. We install version 3.9 because of a bug
    with the most recent version.

For the Python pairinteraction library and the Python GUI, you need a Python 3
distribution (we recommend `Miniconda`_ or `Anaconda`_). The following Python 3
packages have to be installed

.. code-block:: none

    pint psutil pyqt numpy scipy

.. _Miniconda: https://conda.io/miniconda.html
.. _Anaconda: https://www.anaconda.com/distribution/

Build Instructions
""""""""""""""""""

To build for OS X checkout the latest version of pairinteraction
using ``git``

.. code-block:: bash

    $ git clone --single-branch --recursive https://github.com/pairinteraction/pairinteraction.git

.. note::
    Don't forget the ``--recursive`` switch.  Otherwise all the
    submodules will be missing and you won't be able to build
    pairinteraction successfully.  The ``--single-branch`` flag is not
    essential but will speed up the download significantly by omitting
    all other branch except ``master``.

Given that the package ``llvm@3.9`` has been installed via ``homebrew``, we force CMake
to use the Clang C++ compiler by executing the bash commands

.. code-block:: bash

    $ export CXX=/usr/local/opt/llvm@3.9/bin/clang++
    $ export LDFLAGS="-L/usr/local/opt/llvm@3.9/lib -Wl,-rpath,/usr/local/opt/llvm@3.9/lib,-rpath,${CONDA_PREFIX}/lib"

Then proceed with the usual CMake workflow

.. code-block:: bash

    $ cd pairinteraction
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make -j 8

This will build pairinteraction for real and complex matrices.
Afterwards you can start the pairinteraction GUI from the build directory

.. code-block:: bash

    $ ./pairinteraction

In order to use the pairinteraction Python 3 library,
you have to add the build directory to the Python package search path. The pairinteraction C++
library can be used right away.

Code documentation
^^^^^^^^^^^^^^^^^^

To generate the code documentation Doxygen and Sphinx are needed (in
addition to all other build dependencies).  On Ubuntu or Debian it can
be obtained using

.. code-block:: bash

    sudo apt-get install doxygen graphviz python3-sphinx python3-numpydoc

Then checkout the latest version of pairinteraction using `git`

.. code-block:: bash

    $ git clone --single-branch --recursive https://github.com/pairinteraction/pairinteraction.git

.. note::
    Don't forget the ``--recursive`` switch.  Otherwise all the
    submodules will be missing and you won't be able to build
    pairinteraction successfully.  The ``--single-branch`` flag is not
    essential but will speed up the download significantly by omitting
    all other branch except ``master``.

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
