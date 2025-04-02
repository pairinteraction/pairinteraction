# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

# Function to add a test executable
function(add_test_executable TARGET_NAME SOURCE_FILE NPROC)
  # Create the test executable
  add_executable(${TARGET_NAME} ${SOURCE_FILE})
  target_link_libraries(${TARGET_NAME} PRIVATE pairinteraction)

  # Windows-specific post-build command to copy runtime DLLs
  if(WIN32)
    add_custom_command(
      TARGET ${TARGET_NAME}
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:${TARGET_NAME}>
              $<TARGET_FILE_DIR:${TARGET_NAME}>
      COMMAND_EXPAND_LISTS)
  endif()

  # Add tests building and calling the executable
  add_test(NAME gen_${TARGET_NAME} COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config "$<CONFIG>"
                                           --target ${TARGET_NAME} -j${NPROC})
  add_test(NAME ${TARGET_NAME} COMMAND ${TARGET_NAME} --database-dir "${CMAKE_SOURCE_DIR}/data/database" --data-dir
                                       "${CMAKE_SOURCE_DIR}/data")
  set_tests_properties(gen_${TARGET_NAME} PROPERTIES FIXTURES_SETUP f_gen_${TARGET_NAME})
  set_tests_properties(${TARGET_NAME} PROPERTIES FIXTURES_REQUIRED f_gen_${TARGET_NAME})

  # Propagate the name of the test building the executable to the parent scope
  set(GEN_TEST_NAMES
      "${GEN_TEST_NAMES};gen_${TARGET_NAME}"
      PARENT_SCOPE)
endfunction()

# Determine the number of processor cores
include(ProcessorCount)
ProcessorCount(NPROC_DEFAULT)
set(NPROC
    ${NPROC_DEFAULT}
    CACHE STRING "Number of processor cores to use for building the tests")

# Add test executables
add_test_executable(unit_tests unit_tests.cpp ${NPROC})
add_test_executable(test_dipole_operator test_dipole_operator.cpp ${NPROC})
add_test_executable(test_system_atom test_system_atom.cpp ${NPROC})
add_test_executable(test_starkmap test_starkmap.cpp ${NPROC})
add_test_executable(test_pair_potential test_pair_potential.cpp ${NPROC})
add_test_executable(test_ssl test_ssl.cpp ${NPROC})

# Ensure that not more than one test is build at the same time
set_tests_properties(${GEN_TEST_NAMES} PROPERTIES RESOURCE_LOCK compiler_access)
