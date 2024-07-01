.. _Installation:

Installation
============

This guide provides installation instructions for the `pairinteraction` software, catering to users of all levels.
It offers three different installation methods:

- `Using Pip`_: Ideal for most users, this method employs `pip` for straightforward installation from the
  Python Package Index, providing access to the software's Python library and graphical interface.

- `Binary Installers`_: If you are mainly interested in the graphical interface, you can download
  an installer for pairinteraction (Windows, OS X) or use the `Flatpak` package manager (GNU/Linux).

- `Building from Source`_: Aimed at developers and experienced users, this approach involves compiling the software from source,
  detailing necessary dependencies, and offering steps for either an automated build with `poetry` / `pip` or a manual build.
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
    `CMake`_ for running the build system (at least CMake 3.21 is required), `poetry`_ for managing the Python dependencies. In addition, we require a package manager for the C++ dependencies and a compiler for building the C++ backend. In the following, we provide recommendations for different operating systems:

    * For **GNU/Linux**, use the distribution's package manager and the gcc or clang compiler.

    * For **OS X**, use `Homebrew`_ and the clang compiler.

    * For **Windows**, use `VCPKG`_ and `Visual Studio`_.

Dependencies of the C++ backend
    * For **GNU/Linux**, complete lists of dependencies can be found in the Dockerfiles that we use for continuous integration. The Dockerfiles are located in the :github:`docker branch <tree/docker/docker>`.

    * For **OS X**, you can obtain the dependencies from the :github:`macos workflow <tree/master/.github/workflows/macos.yml>`.

    * For **Windows**, you can use VCPKG with :github:`our configuration file <tree/master/vcpkg.json>` to install most dependencies. Further dependencies such as `Intel MKL`_ can be found in the :github:`windows workflow <tree/master/.github/workflows/windows.yml>` and :github:`actions folder <tree/master/.github/actions>` of the pairinteraction repository.

Dependencies of the Python library
    All Python dependencies are listed within the :github:`pyproject.toml <tree/master/pyproject.toml>` file. They are installed automatically when you build the Python library using poetry.

.. _git: https://git-scm.com/download/
.. _CMake: https://cmake.org/download/
.. _poetry: https://python-poetry.org/docs/#installing-with-the-official-installer
.. _Homebrew: https://brew.sh/
.. _VCPKG: https://github.com/microsoft/vcpkg?tab=readme-ov-file#quick-start-windows
.. _Visual Studio: https://visualstudio.microsoft.com/downloads/
.. _Intel MKL: https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html

Automatic Build
^^^^^^^^^^^^^^^

.. note::
    If you do not want to modify the source code and just want to use the most recent version of pairinteraction, you can install pairinteraction directly from the :github:`GitHub <>` repository by running
    ``pip install git+https://github.com/pairinteraction/pairinteraction``. Similarly, you can add the most recent version of pairinteraction to a Python project that is managed by poetry by running ``poetry add git+https://github.com/pairinteraction/pairinteraction``.

After cloning the repository and installing the requirements, you can build and install the software into a local virtual Python environment by running the following command within the pairinteraction repository:

.. code-block:: bash

    poetry install

This will call CMake automatically to build the C++ backend, the Python library, and the graphical user interface. The graphical user interface can be started by executing

.. code-block:: bash

    poetry run start_pairinteraction_gui

To use the Python library, you have to run your python code in the virtual environment created by poetry. This can be done by running ``poetry run python your_script.py``.
Alternatively, you can build and install the software system-wide by running ``pip install -e .`` from the root directory of the pairinteraction repository.

Tests of the Python library and graphical user interface can be run by executing

.. code-block:: bash

    poetry run pytest

Manual Build
^^^^^^^^^^^^

.. note::
    Advanced examples for the usage of CMake to build the software for various operating systems can be found in the :github:`workflows <tree/master/.github/workflows>` directory of the pairinteraction repository.

If you want to build, e.g., the documentation of pairinteraction or have more control over the build process, you can run the tasks that have been executed by poetry manually. For this, you have to first install the Python dependencies manually:

.. code-block:: bash

    poetry export -f requirements.txt > requirements.txt
    pip install -r requirements.txt

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
    cmmake --build . --config RelWithDebInfo

This creates the C++ backend, the Python library, and the graphical user interface. The graphical user interface can be started by executing

.. code-block:: bash

    ./start_pairinteraction_gui

in the build directory.
To use the Python library, you have to extend the Python package search path to accommodate pairinteraction by adding your build directory to ``PYTHONPATH``.
This can be done e.g. by adding the following lines to the top of a Python script:

.. code-block:: python

    import sys
    sys.path.append("/your/path/to/pairinteraction/build")

Running the different build commands manually has the advantage that you can pass additional options to the build system. For example, you can disable the graphical user interface by running CMake with ``cmake -DWITH_GUI=OFF ..`` (the general format to set an option is ``-D<OPTION_NAME>=<VALUE>``).
A full list of build options is provided in the following:

+---------------------+--------------------------------------+---------+
| Option              | Effect                               | Default |
+=====================+======================================+=========+
| ``WITH_BACKEND``    | Build with C++ backend               | ON      |
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
| ``WITH_DMG``        | Generate a binary DMG file (Mac OS X | OFF     |
|                     | only)                                |         |
+---------------------+--------------------------------------+---------+
| ``WITH_COVERAGE``   | Generate code coverage report        | OFF     |
+---------------------+--------------------------------------+---------+
| ``WITH_LTO``        | Build with link-time optimization    | OFF     |
+---------------------+--------------------------------------+---------+
| ``WITH_CLANG_TIDY`` | Run the C++ linter tool Clang-Tidy   | OFF     |
|                     | during compilation                   |         |
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
| ``livehtml`` | Build the Sphinx documentation and        | ``WITH_DOC=ON``      |
|              | show it in a web browser                  |                      |
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

    cmake --build . --config Release --target pairinteraction_backend_test
    ctest -V -C Release -R pairinteraction_backend_test

However, before pushing your changes, you should always run the full test suite to ensure that your changes do not break other parts of the software. The ``--config Release`` and ``-C Release`` options tell the tools to build and test the release version of the software if a multi-configuration generator is used. For further explanations on the build type, see the next section.

**5. Debugging with GDB**

For tracking down errors like segmentation faults, running a debug build with the GNU Debugger `GDB` can be very helpful.

If CMake uses a multi-configuration generator (e.g., Ninja Multi-Config, Visual Studio Generators), you can build the software with debug symbols by using the ``--config Debug`` option. Afterwards, you can execute the build with GDB. For example:

.. code-block:: bash

    cmake -G"Ninja Multi-Config" -DCMAKE_CXX_FLAGS="-fuse-ld=mold" ..
    cmake --build . --config Debug --target pairinteraction_backend_test
    gdb -ex r --args pairinteraction_backend/Debug/pairinteraction_backend_test

If you are using a single-configuration generator (e.g., Unix Makefiles), you must specify the build type directly:

.. code-block:: bash

    cmake -DCMAKE_BUILD_TYPE=Debug ..
    cmake --build . --target pairinteraction_backend_test
    gdb -ex r --args pairinteraction_backend/pairinteraction_backend_test

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
