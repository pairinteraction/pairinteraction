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