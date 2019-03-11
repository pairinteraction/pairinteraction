# Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
#
# This file is part of the pairinteraction library.
#
# The pairinteraction library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The pairinteraction library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.

# Find Numpy
#
# NUMPY_INCLUDE_DIR   - directory with all header files

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include (FindPackageHandleStandardArgs)

if(NOT PYTHON_EXECUTABLE)
  find_package(PythonInterp 3)
endif()

if(PYTHON_EXECUTABLE)
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
    "from numpy import get_include; print(get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
    "from numpy.version import version; print(version)"
    OUTPUT_VARIABLE NUMPY_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  message(STATUS "Python not found, not finding Numpy")
endif()

find_path( NUMPY_INCLUDE_DIR
  NAMES numpy/arrayobject.h
  HINTS ${NUMPY_INCLUDE_PATH} ${PYTHON_INCLUDE_PATH} )

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args( Numpy DEFAULT_MSG NUMPY_INCLUDE_DIR )

mark_as_advanced( NUMPY_INCLUDE_DIR )
