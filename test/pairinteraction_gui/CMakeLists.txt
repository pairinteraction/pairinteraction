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
  add_test(NAME "python.${UT_TARGET}" COMMAND ${Python3_EXECUTABLE} ${UT_SOURCE})
  set_tests_properties("python.${UT_TARGET}" PROPERTIES
    ENVIRONMENT "${UT_ENVIRONMENT};PYTHONPATH=$<SHELL_PATH:${CMAKE_BINARY_DIR};${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE};$ENV{PYTHONPATH}>")
endfunction(PYTHON_TEST)

# Add Python test cases
if (WITH_PYTHON)
  # if(WITH_GUI AND NOT (APPLE AND DEFINED ENV{CI}))
  #   python_test(TARGET gui SOURCE gui.py)
  #   file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/reference_data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  # endif() # TODO somehow this test times out on my machine
endif()
