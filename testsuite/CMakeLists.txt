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

# Build tests only when testing
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL ON)

find_package(PythonInterp 3 REQUIRED)

# To add test cases easier we use the following function
function(PYTHON_TEST)
  # Parse arguments
  cmake_parse_arguments(UT "" "TARGET" "SOURCE" ${ARGN} )
  # Copy test to binary directory
  configure_file(${UT_SOURCE} ${CMAKE_CURRENT_BINARY_DIR}/${UT_SOURCE})
  # Tests are executed in the root directory
  add_test(NAME "python.${UT_TARGET}" COMMAND ${PYTHON_EXECUTABLE} ${UT_SOURCE})
  set_tests_properties("python.${UT_TARGET}"
    PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR};${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}")
endfunction(PYTHON_TEST)

# Build library before testing
add_dependencies(check _pireal _picomplex)

# Add Python test cases
python_test(TARGET rotation SOURCE rotation.py)
python_test(TARGET cache SOURCE cache.py)
python_test(TARGET perturbation SOURCE perturbation.py)
python_test(TARGET states SOURCE state.py)
python_test(TARGET compute SOURCE compute.py)
python_test(TARGET quantum_defect SOURCE quantum_defect.py)
python_test(TARGET wavefunction SOURCE wavefunction.py)
python_test(TARGET symmetries SOURCE symmetries.py)
python_test(TARGET integration SOURCE integration.py)
if(NOT MSVC AND NOT (APPLE AND DEFINED ENV{TRAVIS}) AND NOT WITH_CLANG_TIDY) # timeout
    python_test(TARGET parallelization SOURCE parallelization.py)
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/integration_test_referencedata.pickle
  ${CMAKE_CURRENT_BINARY_DIR}/integration_test_referencedata.pickle
  COPYONLY)
