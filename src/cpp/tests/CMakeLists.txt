# SPDX-FileCopyrightText: 2024 PairInteraction Developers
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

# Get the path to the certifi CA bundle using Python
execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import certifi; print(certifi.where(), end='')"
    OUTPUT_VARIABLE CERTIFI_CA_BUNDLE_PATH
    RESULT_VARIABLE CERTIFI_RESULT
    ERROR_VARIABLE CERTIFI_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
)
if(NOT CERTIFI_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to query certifi.where() using ${Python_EXECUTABLE}: ${CERTIFI_ERROR}")
endif()
if(NOT EXISTS "${CERTIFI_CA_BUNDLE_PATH}")
    message(FATAL_ERROR "certifi CA bundle does not exist at ${CERTIFI_CA_BUNDLE_PATH}")
endif()

# Add test executables
add_test_executable(unit_tests unit_tests.cpp ${NPROC})
add_test_executable(test_dipole_operator test_dipole_operator.cpp ${NPROC})
add_test_executable(test_system_atom test_system_atom.cpp ${NPROC})
add_test_executable(test_starkmap test_starkmap.cpp ${NPROC})
add_test_executable(test_pair_potential test_pair_potential.cpp ${NPROC})
add_test_executable(test_ssl test_ssl.cpp ${NPROC})
target_compile_definitions(test_ssl PRIVATE
    CERTIFI_CA_BUNDLE_PATH="${CERTIFI_CA_BUNDLE_PATH}"
)

# Ensure that not more than one test is build at the same time
set_tests_properties(${GEN_TEST_NAMES} PROPERTIES RESOURCE_LOCK compiler_access)
