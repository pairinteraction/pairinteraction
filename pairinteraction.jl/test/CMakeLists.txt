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

# Build tests only when testing
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL ON)

function(JULIA_TEST)
  # Parse arguments
  cmake_parse_arguments(UT "" "TARGET" "SOURCE" ${ARGN} )
  # Copy test to binary directory
  configure_file(${CMAKE_SOURCE_DIR}/${LIBNAME}.jl/test/${UT_SOURCE} ${CMAKE_CURRENT_BINARY_DIR}/${UT_SOURCE})
  # Tests are executed in the root directory
  add_test(NAME "julia.${UT_TARGET}" COMMAND ${Julia_EXECUTABLE} ${UT_SOURCE})
  set_tests_properties("julia.${UT_TARGET}"
    PROPERTIES ENVIRONMENT "JULIA_LOAD_PATH=${CMAKE_BINARY_DIR}/${LIBNAME}.jl/src/:")
endfunction()

if (WITH_JULIA)
  add_dependencies(check pireal_jl picomplex_jl)
  find_package(JuliaInterp REQUIRED)
  julia_test(TARGET field_combination SOURCE test_field_combination.jl)
  julia_test(TARGET green_tensor SOURCE test_green_tensor.jl)
  julia_test(TARGET pair_state SOURCE test_pair_state.jl)
  julia_test(TARGET quantum_defect SOURCE test_quantum_defect.jl)
  julia_test(TARGET state SOURCE test_state.jl)
endif()
