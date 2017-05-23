# Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Find Numpy
#
# NUMPY_INCLUDE_DIR   - directory with all header files

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include (FindPackageHandleStandardArgs)

if(NOT PYTHON_EXECUTABLE)
  find_package(PythonInterp)
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
