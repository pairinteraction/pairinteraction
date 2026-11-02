.. _system_setup:

Setup of the Development Environment
====================================

This section described how to set up your system for a :ref:`manual build <manual>` of the C++ backend from source. This
is highly platform dependent, and is described below for each operating system individually. In order to create a Python
environment for the Python backend, refer to our :ref:`Python <python_setup>` setup instructions.

Windows
-------

In order to be able to compile the source code, you have to install the following tools and dependencies:

Build tools
    - CMake_ for running the build system (at least CMake 3.21 is required)
    - VCPKG_ for managing the C++ packages
    - `Visual Studio`_ or MinGW_ as a compiler. We recommend `Visual Studio`_ and will use it in further descriptions.

You can use VCPKG with :github:`our configuration file <tree/master/vcpkg.json>` to install most C++ dependencies.
Further dependencies such as `Intel oneAPI MKL`_ and `Intel oneAPI TBB`_ can be found in the :github:`github workflow
<tree/master/.github/workflows/cpp-backend.yml>` and :github:`actions folder <tree/master/.github/actions>` of the
PairInteraction repository.

In addition, you need to adjust your path environment variable if you want to use certain tools from the command line.
In order to smoothly run all the commands described on this page, add the following paths to your path environment
variable:

======================= ========================================================
Tool                    Path to add
======================= ========================================================
Cmake                   ``C:\\path\to\cmake\bin``
VCPKG                   ``C:\\path\to\vcpkg``
clang-tidy clang-format ``C:\\path\to\VisualStudio\Community\VC\Tools\Llvm\bin``
======================= ========================================================

You can either adjust your environmental variables by using the GUI provided by the operating system, or immediately
from the command line. If you want to use the Windows GUI, which we highly recommend, enter "environment properties"
into the search bar and hit Enter. In the System Properties window, click "Environment Variables." Then click on
environment variables. If you want to setup the environment for all Users, search in the "System Variables", otherwise
in the "User Variables for UserName" for the variable `Path`. Click on `Path` and then on `Edit`. You now see a list of
all paths that have been added to `path`. Click on `Add` and enter the paths to the tools you want to use in the command
prompt. If you want to only change your environment variabels temporarily, you can do so by typing

.. code-block:: bash

    set PATH=%PATH%;C:\your\path\here

.. warning::

    Even though the command `setx` allows you to permanently override your path variable from the command line, we
    highly disencourage you from using it, as it might permamently delete parts of your system environment variable.

OS X
----

Before compiling the source code, you have to install the following tools and dependencies:

Build tools
    - CMake_ for running the build system (at least CMake 3.21 is required)
    - Homebrew_ and the clang compiler.

You can obtain the C++ dependencies from the :github:`github workflow <tree/master/.github/workflows/cpp-backend.yml>`.

GNU/Linux
---------

Before compiling the source code, you have to install the following tools and dependencies:

Build tools
    - CMake_ for running the build system (at least CMake 3.21 is required)
    - the distribution's package manager for managing C++ packages
    - gcc or clang compiler for compiling C++

A complete lists of all C++ dependencies dependencies can be found in the Dockerfiles that we use for continuous
integration. The Dockerfiles are located in the :github:`docker branch <tree/docker/docker>`.

.. _cmake: https://cmake.org/download/

.. _homebrew: https://brew.sh/

.. _intel oneapi mkl: https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html

.. _intel oneapi tbb: https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb-download.html

.. _intel repository: https://www.intel.com/content/www/us/en/docs/oneapi/installation-guide-linux/2023-0/apt.html

.. _mingw: https://www.mingw-w64.org/downloads/

.. _vcpkg: https://github.com/microsoft/vcpkg?tab=readme-ov-file#quick-start-windows

.. _visual studio: https://visualstudio.microsoft.com/downloads/

Tips independent of the OS
--------------------------

.. note::

    `Intel oneAPI MKL`_ is an optional dependency that provides optimized linear algebra routines and the FEAST
    eigensolver. If this dependency is available, it is important that a compatible version of `Intel oneAPI TBB`_ is
    available as well. For example, on Debian, the package ``intel-oneapi-mkl-devel-2023.1.0`` and
    ``intel-oneapi-tbb-devel-2021.13`` are compatible that can be installed using APT from the `Intel repository`_. To
    allow the build system to find these oneAPI libraries, one has to set the ``CMAKE_PREFIX_PATH`` and
    ``LD_LIBRARY_PATH`` environment variables. To do so, the libraries provide scripts that can be sourced before
    running CMake. On Debian, ``source /opt/intel/oneapi/mkl/latest/env/vars.sh`` and ``source
    /opt/intel/oneapi/tbb/latest/env/vars.sh`` will set the environment variables.

.. note::

    Advanced examples for the usage of CMake to build the software for various operating systems can be found in the
    :github:`workflows <tree/master/.github/workflows>` directory of the PairInteraction repository.
