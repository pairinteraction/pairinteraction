# pairinteraction

The pairinteraction software calculates pair potential for Rydberg
atoms in external fields.

# Dependencies

## calc

- C++11
- Eigen3
- MPI (OpenMPI or MSMPI)
- Boost (`filesystem` and `system`)
- sqlite3
- jsoncpp

## gui

- Python3
- psutil
- PyQt5
- Pint
- numpy
- scipy

# Windows

## Install

So far we tackled cross-compilation from OpenSUSE 13.2 to Windows
32-bit and 64-bit.  To install the MinGW toolchain install the
following packages:

    mingw32-boost-devel
    mingw32-cross-binutils
    mingw32-cross-breakpad-tools
    mingw32-cross-cpp
    mingw32-cross-gcc
    mingw32-cross-gcc-c++
    mingw32-cross-gcc-fortran
    mingw32-filesystem
    mingw32-headers
    mingw32-libboost_filesystem
    mingw32-libboost_system
    mingw32-libgcc_s_sjlj1
    mingw32-libgfortran3
    mingw32-libquadmath0
    mingw32-libsqlite3-0
    mingw32-libstdc++6
    mingw32-libwinpthread1
    mingw32-runtime
    mingw32-sqlite-devel
    mingw32-winpthreads-devel

    mingw64-boost-devel
    mingw64-cross-binutils
    mingw64-cross-breakpad-tools
    mingw64-cross-cpp
    mingw64-cross-gcc
    mingw64-cross-gcc-c++
    mingw64-cross-gcc-fortran
    mingw64-filesystem
    mingw64-headers
    mingw64-libboost_filesystem
    mingw64-libboost_system
    mingw64-libgcc_s_seh1
    mingw64-libgfortran3
    mingw64-libquadmath0
    mingw64-libsqlite3-0
    mingw64-libstdc++6
    mingw64-libwinpthread1
    mingw64-runtime
    mingw64-sqlite-devel
    mingw64-winpthreads-devel

This is not the whole story; furthermore we need
[Eigen3](http://eigen.tuxfamily.org/) and
[JsonCpp](https://github.com/open-source-parsers/jsoncpp), both of
which we download the latest release version and compile from source.

For Eigen3, we have to comment out two lines in the file
`cmake/EigenDetermineOSVersion.cmake`.  Otherwise CMake refuses to
generate Makefiles.

    list(GET ver_list 0 _major)
    list(GET ver_list 1 _minor)

For JsonCpp we have to disable the testsuite, by setting the option
`JSONCPP_WITH_TESTS` to `OFF`.  It needs extra libraries we don't want
to install and in addition to that, I believe that the tests pass
anyway.

---

For the GUI you need to install a Python distribution with the
dependencies listed above.  We have tested the installation using the
Anaconda and Miniconda distributions; see instructions below.

### Anaconda

Download [Anaconda3](https://www.continuum.io/downloads)

Install dependencies:

    pip install pint pyqt5

### Miniconda

Download [Miniconda3](http://conda.pydata.org/miniconda.html)

Install dependencies:

    conda install numpy scipy
    pip install psutil pint pyqt5

## Known Bugs

The program crashes on 32-bit Windows when the basis is chosen too
large.  This is due to the fact that 32-bit Windows limits the user
address space for an individual process to 2GB.  If you run
`pairinteraction` from CMD you will see that the `std::bad_alloc`
exception was thrown.  (This could probably be resolved by compiling
for 64-bit, but our current cross-compilation toolchain does not
support it.)