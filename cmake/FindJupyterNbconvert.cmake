# Copyright (c) 2019 Sebastian Weber, Henri Menke. All rights reserved.
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

# Find jupyter-nbconvert
#
# JUPYTER_NBCONVERT_EXECUTABLE

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include(FindPackageHandleStandardArgs)

if(NOT PYTHON_EXECUTABLE)
  find_package(PythonInterp 3)
endif()

if(PYTHON_EXECUTABLE)
  # find ipython3
  find_program(JUPYTER_NBCONVERT_EXECUTABLE NAMES jupyter-nbconvert)

  # Get the version
  execute_process(
    COMMAND "${JUPYTER_NBCONVERT_EXECUTABLE}" --version
    OUTPUT_VARIABLE JUPYTER_NBCONVERT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
else()
  message(STATUS "Python not found, not finding IPython")
endif()

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args(JupyterNbconvert
  REQUIRED_VARS JUPYTER_NBCONVERT_EXECUTABLE
  VERSION_VAR JUPYTER_NBCONVERT_VERSION
  )

mark_as_advanced(JUPYTER_NBCONVERT_EXECUTABLE)
