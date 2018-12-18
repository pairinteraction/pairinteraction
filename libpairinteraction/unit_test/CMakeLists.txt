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

# Add current directory to search path
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Find the Boost unit test framework
find_package(Boost REQUIRED COMPONENTS unit_test_framework)
include_directories(${Boost_INCLUDE_DIRS})
list(APPEND LIBRARIES ${Boost_LIBRARIES})

# To add test cases easier we use the following function
function(UNIT_TEST)
  # Parse arguments
  cmake_parse_arguments(UT "" "TARGET" "SOURCE" ${ARGN} )
  # Add target for test and link with libraries
  add_executable(${UT_TARGET} ${UT_SOURCE})
  if(MSVC)
    set_source_files_properties(${UT_SOURCE} PROPERTIES COMPILE_FLAGS "/W3 /DNOMINMAX")
  else()
    set_source_files_properties(${UT_SOURCE} PROPERTIES COMPILE_FLAGS "-Wall -Wextra -pedantic -Werror")
  endif()
  target_link_libraries(${UT_TARGET} ${LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  target_compile_definitions(${UT_TARGET} PRIVATE "BOOST_TEST_DYN_LINK=1")
  # Tests are executed in the root directory
  add_test(NAME ${UT_TARGET} COMMAND ${UT_TARGET}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  add_dependencies(check ${UT_TARGET})
endfunction(UNIT_TEST)

# Add the tests
unit_test(TARGET sqlite_test SOURCE sqlite_test.cpp)
unit_test(TARGET quantum_defect_test
  SOURCE quantum_defect_test.cpp ../QuantumDefect.cpp)
unit_test(TARGET wavefunction_test
  SOURCE wavefunction_test.cpp ../QuantumDefect.cpp ../Wavefunction.cpp)
unit_test(TARGET conf_parser_test
  SOURCE conf_parser_test ../ConfParser.cpp)
unit_test(TARGET integration_test
  SOURCE integration_test ../State.cpp ../QuantumDefect.cpp ../Wavefunction.cpp
  ../WignerD.cpp ../MatrixElementCache.cpp ../SystemOne.cpp ../SystemTwo.cpp ../GreenTensor.cpp)
unit_test(TARGET cache_test SOURCE cache_test.cpp)
unit_test(TARGET utils_test SOURCE utils_test.cpp)

# Copy test dependencies
add_custom_command(TARGET integration_test PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/integration_test_referencedata.txt ${CMAKE_CURRENT_BINARY_DIR}/integration_test_referencedata.txt
)
