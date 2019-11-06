.. _Installation:

Installation
============

Binary Installers
-----------------

All builds come with an easy-to-use graphical user interface for pair potential calculations, taking into
account electric and magnetic fields. In addition, a Python 3 and C++ library is provided which can be used to
write your own code and have more fine-grained control over what pairinteraction does. Quantities as matrix elements,
eigenenergies, or excitation dynamics can be calculated. For usage examples
visit the :ref:`tutorials <Tutorials>` section of the documentation or
download the :github:`Python <tree/master/doc/sphinx/examples_python>`
and :github:`C++ <tree/master/doc/sphinx/examples_cpp>` example programs.

Packages are available for GNU/Linux, Mac OS X, and Windows through
:github:`GitHub Releases <releases>`. For these packages, the graphical user interface works out-of-the-box.
The different packages were built on the following architectures:

- ``deb`` package: Ubuntu 18.04 amd64
- ``rpm`` package: OpenSUSE Leap 15.0 x86_64
- Mac OS X ``dmg``: Mac OS X 10.13 (Xcode 9.4)
- Windows ``exe``: Compiled with Visual Studio 2015

Alternatively, you can install pairinteraction from the `Python Package Index`_ via pip by calling ``pip install pairinteraction`` from the command line.
We recommend this way of installation for using pairinteraction as a Python 3 library. To get the graphical user interface work with such an installation, the Python
library PyQt5 must be installed as well.

If pairinteraction was installed from the command line, the graphical user interface can be started by executing ``start_pairinteraction_gui``.

.. _Python Package Index: https://pypi.org/project/pairinteraction

Python Library
^^^^^^^^^^^^^^

If pairinteraction was installed via pip or using a Linux package manager, all dependencies are installed
automatically and the library can be used for Python 3 right away.

Otherwise, the dependencies of the pairinteraction library must be installed manually. Assuming that the Python 3
distribution `Miniconda`_ or `Anaconda`_ is used, the dependencies can be installed using conda. Note that it is
important that the installed version of the NumPy library is up-to-date.

.. _Miniconda: https://conda.io/miniconda.html
.. _Anaconda: https://www.anaconda.com/distribution/

.. code-block:: bat

    conda install numpy scipy matplotlib

In addition, the path containing the pairinteraction library must be added manually to the Python package search path.
If the pairinteraction software was installed to its default location, this
can be done by adding the following code block to the top of a Python script:

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
    If you built from source, you have to extend the Python package search path to
    accomodate pairinteraction by adding your build directory to ``PYTHONPATH``. This can be done e.g. by
    adding the following lines to the top of a Python script:
    
    .. code-block:: python

        import sys
        sys.path.append("/your/path/to/pairinteraction/build")
   
Requirements
^^^^^^^^^^^^

The following tools and libraries, including header files, are required
for compiling the source code.

Git
    Git is a version control system used to track changes.

CMake 3.9 or later
    The build system is based on CMake.

C++ Backend and Python Interface
""""""""""""""""""""""""""""""""

C++ Compiler
    C++14 capable C++ compiler

Sqlite3
   SQLite is a self-contained, high-reliability, embedded,
   full-featured, public-domain, SQL database engine.

Boost
    A library providing advanced C++ features.

GSL
    The GNU Scientific Library (GSL) is a numerical library.

SWIG 3.0 or later
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

The build system uses CMake and has some configuration switches. These are

+---------------------+--------------------------------------+---------+
| Option              | Effect                               | Default |
+=====================+======================================+=========+
| ``WITH_BACKEND``    | Build with C++ backend               | ON      |
+---------------------+--------------------------------------+---------+
| ``WITH_PYTHON``     | Build with SWIG Python interface     | ON      |
+---------------------+--------------------------------------+---------+
| ``WITH_GUI``        | Build with Python GUI                | ON      |
+---------------------+--------------------------------------+---------+
| ``WITH_DATABASE``   | Store quantum defect database so     | ON      |
|                     | that it can be accessed by the GUI   |         |
+---------------------+--------------------------------------+---------+
| ``WITH_DOCS``       | Generate documentation               | ON      |
+---------------------+--------------------------------------+---------+
| ``WITH_DMG``        | Generate a DMG file (Mac OS X only)  | OFF     |
+---------------------+--------------------------------------+---------+
| ``WITH_COVERAGE``   | Generate code coverage report        | OFF     |
+---------------------+--------------------------------------+---------+
| ``WITH_LTO``        | Build with link-time optimization    | OFF     |
+---------------------+--------------------------------------+---------+
| ``WITH_GSL``        | Use the GNU scientific library for   | ON      |
|                     | Whittaker functions [#]_             |         |
+---------------------+--------------------------------------+---------+
| ``WITH_CLANG_TIDY`` | Run Clang-Tidy during compilation    | OFF     |
+---------------------+--------------------------------------+---------+

.. [#] This mode activates the extension for calculating radial wave
       functions using Whittaker functions. If pairinteraction
       is built in this mode, any derived work has to be licensed under
       GPL v3, because of the GSL being distributed under GPL. 

These options can be passed directly to ``cmake``, i.e.

.. code-block:: bash

    $ cmake -DWITH_GUI=OFF -DWITH_DATABASE=OFF -DWITH_DOCS=OFF ..

This way we can only build the C++ backend with the Python interface.

Ubuntu 18.04
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

    patterns-devel-C-C++-devel_C_C++ sqlite3 sqlite3-devel libboost_filesystem1_66_0-devel libboost_program_options1_66_0-devel libboost_serialization1_66_0-devel libboost_system1_66_0-devel libboost_test1_66_0-devel gsl-devel swig python3 python3-devel python3-numpy python3-numpy-devel python3-scipy

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

    cmake git gsl swig libomp

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

    sudo apt-get install doxygen graphviz python3-sphinx python3-numpydoc jupyter-nbconvert pandoc python3-ipykernel python3-matplotlib

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
