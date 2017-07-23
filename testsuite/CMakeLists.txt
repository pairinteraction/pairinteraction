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
  # Tests are executed in the "calc" directory
  add_test(NAME ${UT_TARGET} COMMAND ${PYTHON_EXECUTABLE} ${UT_SOURCE})
  set_tests_properties(${UT_TARGET}
    PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}")
endfunction(PYTHON_TEST)

# Copy python tests
add_custom_target(copy_python_tests
  COMMAND ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )
add_dependencies(check copy_python_tests)

# Add test cases
python_test(TARGET states SOURCE state.py)