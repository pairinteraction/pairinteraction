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

# Find Julia
#
# Julia_EXECUTABLE - the Julia interpreter

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include (FindPackageHandleStandardArgs)

find_program(Julia_EXECUTABLE NAMES julia)

if(Julia_EXECUTABLE)
  execute_process(
    COMMAND "${Julia_EXECUTABLE}" --startup-file=no --version
    OUTPUT_VARIABLE Julia_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  string(REGEX REPLACE ".*([0-9]+\\.[0-9]+\\.[0-9]+).*" "\\1" Julia_VERSION "${Julia_VERSION}")
endif()

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args(JuliaInterp
    REQUIRED_VARS Julia_EXECUTABLE
    VERSION_VAR Julia_VERSION
)
mark_as_advanced(Julia_EXECUTABLE)
