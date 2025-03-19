.. _manual:

Build using cmake
=================

Build process
-------------

If you want to build only the C++ part and want to have more control over the build process, you can manually run the tasks that have been automatically executed by `pip`_ in the :ref:`automatic build <automatic>`.
For this, you have to first install the python build dependencies for your :ref:`python environment <python_setup>` manually.

If you want to use mkl you should also run ``pip install mkl mkl-devel``.

You can then build the software with standard CMake commands:

.. code-block:: bash

    mkdir build
    cd build
    cmake ..
    cmake --build .

Make sure that on a :ref:`Windows platform <system_setup>`, you have set the correct environmental variables.
If you are using a multi-configuration generator for CMake, you must also specify a configuration at the build command.


Running the different build commands manually has the advantage that you can pass additional options to the build system. For example, you can enable the code coverage by running CMake with ``cmake -DWITH_COVERAGE=ON ..`` (the general format to set an option is ``-D<OPTION_NAME>=<VALUE>``).
A full list of build options is provided in the following:

+---------------------+--------------------------------------+---------+
| Option              | Effect                               | Default |
+=====================+======================================+=========+
| ``WITH_COVERAGE``   | Generate code coverage report [1]_   | OFF     |
+---------------------+--------------------------------------+---------+

.. [1] This mode implies building the debug version of the software.

Moreover, executing the commands manually allows for running additional targets.
For example, you can use the ``doxygen`` target to build the C++ `doxygen documentation <https://www.pairinteraction.org/pairinteraction/doxygen/html/index.html>`_ by executing ``cmake --build . --target doxygen``.
In contrast, if you use `pip`_ to build the software, only the default target for building the library is executed.
In the following, a list of all available targets is provided.
Note that some targets require specific build options to be enabled in addition to the default options, and have varying names depending on the platform.

+------------------+----------------------------+------------------------------------------+----------------------+
| Target (Windows) | Target (OS X and Unix)     |Task                                      | Requirement          |
+==================+============================+==========================================+======================+
| ``ALL_BUILD``    | ``all``                    | Build the software (default target)      |                      |
+------------------+----------------------------+------------------------------------------+----------------------+
| ``RUN_TESTS``    | ``test``                   | Run the C++ tests                        |                      |
|                  |                            | (without any python tests,               |                      |
|                  |                            | use the automatic build above for this)  |                      |
+------------------+----------------------------+------------------------------------------+----------------------+
| ``DOXYGEN``      | ``doxygen``                | Build the Doxygen documentation          |                      |
|                  |                            | in ``src/cpp/docs``                      |                      |
+------------------+----------------------------+------------------------------------------+----------------------+

.. note::
    Since ``Visual Studio 17 2022`` is a multi-configuration generator, you have to specify a configuration, otherwise the build will not work. You can choose between ``Release``, ``RelWithDebInfo``, ``Debug`` and ``MinSizeRel``. They are optimizing different purposes:

    +--------------+-----------------------------------------------+
    |Configuration | Purpose                                       |
    +==============+===============================================+
    |Release       |Performance optimized, no debug information    |
    +--------------+-----------------------------------------------+
    |Debug         |Includes debugging information, no optimization|
    +--------------+-----------------------------------------------+
    |RelWithDebInfo|Combines optimization with Debug information   |
    +--------------+-----------------------------------------------+
    |MinSizeRel    |Optimizes binary size of output                |
    +--------------+-----------------------------------------------+

    When testing, you should use the same configuration as in the build, for example

    .. code-block:: bash

      cmake --build . --config RelWithDebInfo
      ctest -C RelWithDebInfo


In addition, a number of options are typically available for the native build tool that is called by CMake.
For example, you can pass the ``-j num_jobs`` option to the native build tool to enable parallel compilation,
where ``num_jobs`` specifies the maximal number of jobs that will be run. Setting ``num_jobs`` to the number of available
processors can speed up the compilation process significantly.

.. code-block:: bash

    cmake --build . -j 8


Tips and Tricks
---------------

**1. Compiler Optimizations**

To speed up the software, you can pass optimization flags to the compiler by setting the `CXXFLAGS` environment variable before running CMake. For example, the following bash command sets the environment variable under GNU/Linux, enabling several optimizations at once for the gcc compiler:

.. code-block:: bash

    export CXXFLAGS="-march=x86-64-v3"

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

**6. Changing the log level**

We use the `spdlog`_ library for logging. The log level can be set by the environment variable `SPDLOG_LEVEL`. Possible values are `info` (the default), `debug`, `warn`, and `error`.

.. _spdlog: https://github.com/gabime/spdlog/

**7. Debugging with GDB**

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
| Basics                  |                                                                  |
+=========================+==================================================================+
| ``help COMMAND``        | Display help for the given COMMAND                               |
+-------------------------+------------------------------------------------------------------+
| ``q``                   | Quit the debugger                                                |
+-------------------------+------------------------------------------------------------------+

+---------------------------+------------------------------------------------------------------+
| Investigating a backtrace |                                                                  |
+===========================+==================================================================+
| ``bt``                    | Display a backtrace of the call stack                            |
+---------------------------+------------------------------------------------------------------+
| ``frame NUMBER``          | Select the frame with the given NUMBER on the call stack         |
+---------------------------+------------------------------------------------------------------+
| ``up`` / ``down``         | Select one frame up or down from the currently selected frame    |
+---------------------------+------------------------------------------------------------------+
| ``list``                  | Display code around the selected frame                           |
+---------------------------+------------------------------------------------------------------+
| ``p EXPR``                | Display the value of EXPR                                        |
+---------------------------+------------------------------------------------------------------+

+----------------------------------+------------------------------------------------------------------+
| Debugging with multiple threads  |                                                                  |
+==================================+==================================================================+
| ``info threads``                 | Display all threads running in the program, the first            |
|                                  | field is the thread number                                       |
+----------------------------------+------------------------------------------------------------------+
| ``thread NUMBER``                | Select the thread with the given NUMBER                          |
+----------------------------------+------------------------------------------------------------------+

+--------------------------+------------------------------------------------------------------+
| Breakpoints and stepping |                                                                  |
+==========================+==================================================================+
| ``b FUNCTIONNAME``       | Set breakpoint at FUNCTIONNAME                                   |
+--------------------------+------------------------------------------------------------------+
| ``delete FUNCTIONNAME``  | Delete breakpoint at FUNCTIONNAME                                |
+--------------------------+------------------------------------------------------------------+
| ``c``                    | Continue executing the program until the next breakpoint         |
+--------------------------+------------------------------------------------------------------+
| ``n``                    | Execute next source-code line, stepping over function calls      |
+--------------------------+------------------------------------------------------------------+
| ``s``                    | Execute next source-code line, stepping into function calls      |
+--------------------------+------------------------------------------------------------------+

.. _gdb's commands: http://www.unknownroad.com/rtfm/gdbtut/gdbtoc.html
.. _pip: https://pypi.org/project/pip/
