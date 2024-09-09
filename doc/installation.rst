.. _Installation:

Installation
============

This guide provides installation instructions for the `pairinteraction` software, catering to users of all levels.
It offers three different installation methods:

- `Using Pip`_: Ideal for most users, this method employs `pip`_ for straightforward installation from the
  Python Package Index, providing access to the software's Python library and graphical interface.

- `Binary Installers`_: If you are mainly interested in the graphical interface, you can download
  an installer for pairinteraction (Windows, OS X) or use the `Flatpak` package manager (GNU/Linux).

- `Building from Source`_: Aimed at developers and experienced users, this approach involves compiling the software from source,
  detailing necessary dependencies, and offering steps for either an automated build with `pip`_ / `uv`_ or a manual build.
  The latter allows for additional customization and the execution of different targets, such as generating documentation.

All methods install the graphical user interface of pairinteraction. It allows for calculating pair potential, taking into
account electric and magnetic fields. In addition, all methods except for the binary installers, provide a Python library which can be used to
write your own code and have more fine-grained control over what pairinteraction does. For usage examples
visit the :ref:`tutorials <Tutorials>` section of the documentation.

Using Pip
---------

We recommend installing pairinteraction from the `Python Package Index`_ by calling

.. code-block:: bash

    pip install pairinteraction

from the command line. This gives you access to the pairinteraction Python library and the graphical user interface.
The graphical user interface can be started by executing

.. code-block:: bash

    start_pairinteraction_gui

Binary Installers
-----------------

Alternatively, if you are mainly interested in the graphical user interface, you can download an installer for pairinteraction from :github:`GitHub Releases <releases>` (Windows, OS X) or
use the `Flatpak`_ package manager (GNU/Linux). For the installation of the Flatpak package, you have to `install Flatpak`_ first and
then run ``flatpak install org.pairinteraction.Pairinteraction`` from the command line.

.. _Python Package Index: https://pypi.org/project/pairinteraction
.. _Flatpak: https://flathub.org/apps/org.pairinteraction.Pairinteraction
.. _install Flatpak: https://flathub.org/setup


Building from Source
--------------------

Advanced users, especially those who want to contribute to the development of pairinteraction, can build the software from source. The source code is available on
:github:`GitHub <>` and can be cloned from the pairinteraction repository using the following `git`_ command:

.. code-block:: bash

    git clone --single-branch https://github.com/pairinteraction/pairinteraction.git

The ``--single-branch`` flag is not essential but will speed up the download significantly by omitting all other branches except master.

Requirements
^^^^^^^^^^^^

Before compiling the source code, you have to install the following tools and dependencies:

Build tools
    `CMake`_ for running the build system (at least CMake 3.21 is required), `pip`_ / `uv`_ for managing the Python dependencies. In addition, we require a package manager for the C++ dependencies and a compiler for building the C++ backend. In the following, we provide recommendations for different operating systems:

    * For **GNU/Linux**, use the distribution's package manager and the gcc or clang compiler.

    * For **OS X**, use `Homebrew`_ and the clang compiler.

    * For **Windows**, use `VCPKG`_ and `Visual Studio`_.

Dependencies of the C++ backend
    * For **GNU/Linux**, complete lists of dependencies can be found in the Dockerfiles that we use for continuous integration. The Dockerfiles are located in the :github:`docker branch <tree/docker/docker>`.

    * For **OS X**, you can obtain the dependencies from the :github:`macos workflow <tree/master/.github/workflows/macos.yml>`.

    * For **Windows**, you can use VCPKG with :github:`our configuration file <tree/master/vcpkg.json>` to install most dependencies. Further dependencies such as `Intel oneAPI MKL`_ and `Intel oneAPI TBB`_ can be found in the :github:`windows workflow <tree/master/.github/workflows/windows.yml>` and :github:`actions folder <tree/master/.github/actions>` of the pairinteraction repository.

    .. note::
        `Intel oneAPI MKL`_ is an optional dependency that provides optimized linear algebra routines and the FEAST eigensolver. If this dependency is available, it is important that a compatible version of `Intel oneAPI TBB`_ is available as well. For example, on Debian, the package ``intel-oneapi-mkl-devel-2023.1.0`` and ``intel-oneapi-tbb-devel-2021.13`` are compatible that can be installed using APT from the `Intel repository`_. To allow the build system to find these oneAPI libraries, one has to set the ``CMAKE_PREFIX_PATH`` and ``LD_LIBRARY_PATH`` environment variables. To do so, the libraries provide scripts that can be sourced before running CMake. On Debian, ``source /opt/intel/oneapi/mkl/latest/env/vars.sh`` and ``source /opt/intel/oneapi/tbb/latest/env/vars.sh`` will set the environment variables.

Dependencies of the Python library
    All Python dependencies are listed within the :github:`pyproject.toml <tree/master/pyproject.toml>` file. They are installed automatically when you build the Python library using `pip`_ / `uv`_.

.. _git: https://git-scm.com/download/
.. _CMake: https://cmake.org/download/
.. _uv: https://pypi.org/project/uv/
.. _pip: https://pypi.org/project/pip/
.. _Homebrew: https://brew.sh/
.. _VCPKG: https://github.com/microsoft/vcpkg?tab=readme-ov-file#quick-start-windows
.. _Visual Studio: https://visualstudio.microsoft.com/downloads/
.. _Intel oneAPI MKL: https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html
.. _Intel oneAPI TBB: https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb-download.html
.. _Intel repository: https://www.intel.com/content/www/us/en/docs/oneapi/installation-guide-linux/2023-0/apt.html

Automatic Build
^^^^^^^^^^^^^^^

**1. Setup**
.. note::
    If you do not want to modify the source code and just want to use the most recent version of pairinteraction, you can install pairinteraction directly from the :github:`GitHub <>` repository by running
    ``pip install git+https://github.com/pairinteraction/pairinteraction``.

We strongly recommend using the `uv`_ tool, to create a virtual environment and build the software inside this virtual environment.
Creating and activating a virtual environment with `uv`_ can be done by running the following commands:

.. code-block:: bash

    uv venv --python=3.8 .venv
    source .venv/bin/activate

In the following, we will describe how to build the software inside this environment using `uv`_.
If not stated otherwise, all commands should be executed from inside the virtual environment and the root directory of the pairinteraction repository.
However, you can also use a different virtual environment manager like `venv` or `conda`, or build the software into your system-wide Python environment by replacing in all following commands ``uv pip`` with ``pip``.

**2. Basic installation**
After cloning the repository and creating the virtual environment, you can build and install the software by running:

.. code-block:: bash

    uv pip install -e .[gui]

This will call CMake automatically to build the C++ backend, the Python library, and the graphical user interface.
By omitting the ``[gui]`` option, you can build the software without installing the additional dependencies needed for the graphical user interface.
The option ``-e`` installs the software in editable mode, which means that changes to the python source code are directly reflected in the installed package.
If you don't want this, you can omit the ``-e`` option to install the current version of the software into the virtual environment.

The graphical user interface can now be started by executing:

.. code-block:: bash

    start_pairinteraction_gui

To use the Python library within your code, you can simply run your python code from inside the virtual environment.

**3. Testing**
Tests of the Python library and graphical user interface can be run by executing:

.. code-block:: bash

    pytest

**4. Build Documentation**
todo

**5. Advanced installation**
Advanced options for developers when building the package:

.. code-block:: bash

    uv pip install --no-build-isolation -Cbuild-dir=build_pip -v -e .

| ``--no-build-isolation``: Avoid re-creations of virtual environments for building the package (to use this you first have to install all build dependencies like ``uv pip install scikit-build-core nanobind numpy``).
| ``-Cbuild-dir=build``: Specify a build directory and reuse it for faster future builds.
| ``-v``: Make the output more verbose.
| ``-e``: Install the package in editable mode (i.e. changes to the python files inside pairinteraction/ are immediately effective).

To install all dependencies without building the package, you can run:
.. code-block:: bash

    uv pip compile pyproject.toml --all-extras > requirements.txt
    uv pip install -r requirements.txt


Manual Build
^^^^^^^^^^^^

.. note::
    Advanced examples for the usage of CMake to build the software for various operating systems can be found in the :github:`workflows <tree/master/.github/workflows>` directory of the pairinteraction repository.

If you want to build only the C++ part and want to have more control over the build process, you can run the tasks that have been executed by `pip`_ / `uv`_ manually.
For this, you have to first install the Python dependencies manually.

Again, we strongly recommend installing the dependencies into a virtual environment using `uv`_:

.. code-block:: bash

    uv venv --python=3.8 .venv
    source .venv/bin/activate

    uv pip compile pyproject.toml --all-extras > requirements.txt
    uv pip install -r requirements.txt

For **GNU/Linux and OS X**, you can then build the software with standard CMake commands:

.. code-block:: bash

    mkdir build
    cd build
    cmake ..
    cmake --build .

For **Windows**, you must specify a visual studio generator, provide a path to the VCPKG toolchain file, and define the build type manually. For example, if you are using Visual Studio 2022, you can build the software with the following commands:

.. code-block:: bash

    mkdir build
    cd build
    cmake -G "Visual Studio 17 2022" -DCMAKE_TOOLCHAIN_FILE=C:/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake ..
    cmake --build . --config RelWithDebInfo

This creates the C++ backend.

Running the different build commands manually has the advantage that you can pass additional options to the build system. For example, you can enable the code coverage by running CMake with ``cmake -DWITH_COVERAGE=ON ..`` (the general format to set an option is ``-D<OPTION_NAME>=<VALUE>``).
A full list of build options is provided in the following:

+---------------------+--------------------------------------+---------+
| Option              | Effect                               | Default |
+=====================+======================================+=========+
| ``WITH_DOXYGEN``    | Generate C++ doxygen documentation   | OFF     |
+---------------------+--------------------------------------+---------+
| ``WITH_DMG``        | Generate a binary DMG file (Mac OS X | OFF     |
|                     | only)                                |         |
+---------------------+--------------------------------------+---------+
| ``WITH_COVERAGE``   | Generate code coverage report [#]_   | OFF     |
+---------------------+--------------------------------------+---------+

.. [#] This mode implies building the debug version of the software.

Moreover, executing the commands manually allows for running additional targets.
For example, you can use the ``doxygen`` target to build the C++ doxygen documentation by executing ``cmake --build . --target doxygen``.
In contrast, if you use `pip`_ / `uv`_ to build the software, only the default target for building the library is executed.
In the following, a list of all available targets is provided.
Note that some targets require specific build options to be enabled in addition to the default options.

+--------------+-------------------------------------------+----------------------+
| Target       | Task                                      | Requirement          |
+==============+===========================================+======================+
| ``all``      | Build the software (default target)       |                      |
+--------------+-------------------------------------------+----------------------+
| ``test``     | Run the C++ tests                         |                      |
|              | (without any python tests,                |                      |
|              | use the automatic build above for this)   |                      |
+--------------+-------------------------------------------+----------------------+
| ``doxygen``  | Build the Doxygen documentation           | ``WITH_DOXYGEN=ON``  |
|              | in ``src/cpp/doc``                        |                      |
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

    cmake --build . -- -j 8

Tips and Tricks
^^^^^^^^^^^^^^^

**1. Compiler Optimizations**

To speed up the software, you can pass optimization flags to the compiler by setting the `CXXFLAGS` environment variable before running CMake. For example, the following bash command sets the environment variable under GNU/Linux, enabling several optimizations at once for the gcc compiler:

.. code-block:: bash

    export CXXFLAGS="-O3 -march=broadwell"

If you are using Windows with Visual Studio, reasonable optimization flags can be set by running the following command in the PowerShell:

.. code-block:: bash

    $env:CXXFLAGS="/Ox /arch:AVX2"

**2. Using a Faster Build System**

Under GNU/Linux, you can use the `ninja` build system and the `mold` linker to reduce the build time by a factor of about 1.5. These tools are typically available in the package repositories of your distribution. For example, on Ubuntu, you can install them by running:

.. code-block:: bash

    sudo apt install ninja-build mold

Then, you can tell CMake to build the software with these tools by running the following commands within the build directory. Note that ninja uses all available processors by default.

.. code-block:: bash

    cmake -G"Ninja Multi-Config" -DCMAKE_CXX_FLAGS="-fuse-ld=mold" ..
    cmake --build .

**3. Using Compiler Caching**

If you delete the build directory because you want to compile a different branch of pairinteraction or use different build options, the compilation has to start from scratch - as long as you do not use a compiler cache like `ccache`. Using this tool has the additional advantage that adding comments to the source code does not trigger a recompilation. It can be installed on many operating systems, e.g., on Ubuntu by running:

.. code-block:: bash

    sudo apt install ccache

To use the tool with CMake, pass ``-DCMAKE_CXX_COMPILER_LAUNCHER=ccache`` to the ``cmake`` command.

**4. Building and Testing Only Parts of the Software**

If you're developing and making changes to specific parts of the software, you can save time by using specific targets to build and test only those parts. You can read off the names of relevant targets from the ``CMakeLists.txt`` files located in the directories where you perform the changes. For example, you can build and test only the C++ backend by running the following commands within the build directory:

.. code-block:: bash

    cmake --build . --config Release --target unit_tests
    ctest -V -C Release -R unit_tests

However, before pushing your changes, you should always run the full test suite to ensure that your changes do not break other parts of the software. The ``--config Release`` and ``-C Release`` options tell the tools to build and test the release version of the software if a multi-configuration generator is used. For further explanations on the build type, see the next section.

**5. Improve the Code Quality with Clang-Tidy and Include-What-You-Use**

Our continues integration system uses the C++ linter tool `clang-tidy` to check the code quality of pull requests and find programming errors. If you have the clang compiler installed, you can run it by yourself during compilation by building the software with the following commands:

.. code-block:: bash

    cmake -DCMAKE_CXX_COMPILER="clang++" -DCMAKE_CXX_CLANG_TIDY="clang-tidy" ..
    cmake --build .

In addition, it is recommended to use the `include-what-you-use` tool to find unnecessary includes in your code. While the tool is not perfect, its suggestions can help to reduce the compilation time. If the tool is installed on your system, you can run it during compilation by executing the following commands:

.. code-block:: bash

    cmake -DCMAKE_CXX_COMPILER="clang++" -DCMAKE_CXX_INCLUDE_WHAT_YOU_USE="iwyu" ..
    cmake --build .

**6. Debugging with GDB**

For tracking down errors like segmentation faults, running a debug build with the GNU Debugger `GDB` can be very helpful.

If CMake uses a multi-configuration generator (e.g., Ninja Multi-Config, Visual Studio Generators), you can build the software with debug symbols by using the ``--config Debug`` option. Afterwards, you can execute the build with GDB. For example:

.. code-block:: bash

    cmake -G"Ninja Multi-Config" -DCMAKE_CXX_FLAGS="-fuse-ld=mold" ..
    cmake --build . --config Debug --target unit_tests
    gdb -ex r --args src/cpp/tests/Debug/unit_tests

If you are using a single-configuration generator (e.g., Unix Makefiles), you must specify the build type directly:

.. code-block:: bash

    cmake -DCMAKE_BUILD_TYPE=Debug ..
    cmake --build . --target unit_tests
    gdb -ex r --args src/cpp/tests/unit_tests

If you have executed a build without GDB, a crash occurred, and a core dump was created, you can load the core dump into GDB:

.. code-block:: bash

    gdb path/to/my/executable path/to/core

After starting the debugger, you can use `GDB's commands`_ to analyze the crash. Some of the most important commands are listed in the tables below.

+-------------------------+------------------------------------------------------------------+
| Basics                                                                                     |
+=========================+==================================================================+
| ``help COMMAND``        | Display help for the given COMMAND                               |
+-------------------------+------------------------------------------------------------------+
| ``q``                   | Quit the debugger                                                |
+-------------------------+------------------------------------------------------------------+

+-------------------------+------------------------------------------------------------------+
| Investigating a backtrace                                                                  |
+=========================+==================================================================+
| ``bt``                  | Display a backtrace of the call stack                            |
+-------------------------+------------------------------------------------------------------+
| ``frame NUMBER``        | Select the frame with the given NUMBER on the call stack         |
+-------------------------+------------------------------------------------------------------+
| ``up`` / ``down``       | Select one frame up or down from the currently selected frame    |
+-------------------------+------------------------------------------------------------------+
| ``list``                | Display code around the selected frame                           |
+-------------------------+------------------------------------------------------------------+
| ``p EXPR``              | Display the value of EXPR                                        |
+-------------------------+------------------------------------------------------------------+

+-------------------------+------------------------------------------------------------------+
| Debugging with multiple threads                                                            |
+=========================+==================================================================+
| ``info threads``        | Display all threads running in the program, the first            |
|                         | field is the thread number                                       |
+-------------------------+------------------------------------------------------------------+
| ``thread NUMBER``       | Select the thread with the given NUMBER                          |
+-------------------------+------------------------------------------------------------------+

+-------------------------+------------------------------------------------------------------+
| Breakpoints and stepping                                                                   |
+=========================+==================================================================+
| ``b FUNCTIONNAME``      | Set breakpoint at FUNCTIONNAME                                   |
+-------------------------+------------------------------------------------------------------+
| ``delete FUNCTIONNAME`` | Delete breakpoint at FUNCTIONNAME                                |
+-------------------------+------------------------------------------------------------------+
| ``c``                   | Continue executing the program until the next breakpoint         |
+-------------------------+------------------------------------------------------------------+
| ``n``                   | Execute next source-code line, stepping over function calls      |
+-------------------------+------------------------------------------------------------------+
| ``s``                   | Execute next source-code line, stepping into function calls      |
+-------------------------+------------------------------------------------------------------+

.. _gdb's commands: http://www.unknownroad.com/rtfm/gdbtut/gdbtoc.html
