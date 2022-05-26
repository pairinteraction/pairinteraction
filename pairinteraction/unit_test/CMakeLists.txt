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

# Build tests only when testing
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL ON)

# Add current directory to search path
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Find the Boost unit test framework
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

# To add test cases easier we use the following function
function(UNIT_TEST)
  # Parse arguments
  cmake_parse_arguments(UT "" "TARGET" "SOURCE" ${ARGN} )
  # Add target for test and link with libraries
  add_executable(${UT_TARGET} ${UT_SOURCE})
  target_compile_features(${UT_TARGET} PRIVATE cxx_std_17)
  set_target_properties(${UT_TARGET} PROPERTIES CXX_EXTENSIONS OFF)
  if(MSVC)
    set_source_files_properties(${UT_SOURCE} PROPERTIES COMPILE_FLAGS "/W3 /DNOMINMAX")
    target_link_libraries(${UT_TARGET} PUBLIC Boost::disable_autolinking Boost::dynamic_linking)
  else()
    set_source_files_properties(${UT_SOURCE} PROPERTIES COMPILE_FLAGS "-Wall -Wextra -pedantic -Werror")
  endif()
  target_link_libraries(${UT_TARGET} PUBLIC pireal Boost::unit_test_framework)
  target_compile_definitions(${UT_TARGET} PRIVATE "BOOST_TEST_DYN_LINK=1")
  # Tests are executed in the root directory
  add_test(NAME "cpp.${UT_TARGET}" COMMAND ${UT_TARGET}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  add_dependencies(check ${UT_TARGET})
endfunction(UNIT_TEST)

# Add the tests
unit_test(TARGET sqlite SOURCE sqlite_test.cpp)
unit_test(TARGET quantum_defect SOURCE quantum_defect_test.cpp)
unit_test(TARGET wavefunction SOURCE wavefunction_test.cpp)
unit_test(TARGET conf_parser SOURCE conf_parser_test.cpp)
unit_test(TARGET integration SOURCE integration_test.cpp)
unit_test(TARGET cache SOURCE cache_test.cpp)
unit_test(TARGET utils SOURCE utils_test.cpp)


# Copy test dependencies
add_custom_command(TARGET integration PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/integration_test_referencedata.txt ${CMAKE_CURRENT_BINARY_DIR}/integration_test_referencedata.txt
)
