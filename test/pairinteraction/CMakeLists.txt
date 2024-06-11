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

# To add test cases easier we use the following function
function(PYTHON_TEST)
  # Parse arguments
  cmake_parse_arguments(UT "" "TARGET" "SOURCE;ENVIRONMENT" ${ARGN} )
  # Copy test to binary directory
  configure_file(${UT_SOURCE} ${CMAKE_CURRENT_BINARY_DIR}/${UT_SOURCE})
  # Tests are executed in the root directory
  add_test(NAME "python.${UT_TARGET}" COMMAND ${Python_EXECUTABLE} ${UT_SOURCE})
  set_tests_properties("python.${UT_TARGET}" PROPERTIES
    ENVIRONMENT "${UT_ENVIRONMENT};PYTHONPATH=$<SHELL_PATH:${CMAKE_BINARY_DIR};${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE};$ENV{PYTHONPATH}>")
endfunction(PYTHON_TEST)

# Add Python test cases
find_package(Python COMPONENTS Interpreter REQUIRED)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
python_test(TARGET model SOURCE test_model.py)
python_test(TARGET simulation SOURCE test_simulation.py)
