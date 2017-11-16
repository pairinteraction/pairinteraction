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

# Find IPython
#
# IPYTHON_EXECUTABLE

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include(FindPackageHandleStandardArgs)

if(NOT PYTHON_EXECUTABLE)
  find_package(PythonInterp 3)
endif()

if(PYTHON_EXECUTABLE)
  # find ipython3
  find_program(IPYTHON_EXECUTABLE NAMES ipython3 ipython)

  # Get the version
  execute_process(
    COMMAND "${IPYTHON_EXECUTABLE}" --version
    OUTPUT_VARIABLE IPYTHON_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
else()
  message(STATUS "Python not found, not finding IPython")
endif()

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args(IPython
  REQUIRED_VARS IPYTHON_EXECUTABLE
  VERSION_VAR IPYTHON_VERSION
  )

mark_as_advanced(IPYTHON_EXECUTABLE)
