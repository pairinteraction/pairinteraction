# Add executable "unit_tests"
add_executable(unit_tests
    unit_tests.cpp
)
target_link_libraries(unit_tests PRIVATE pairinteraction)

if (WIN32)
    add_custom_command(TARGET unit_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:unit_tests> $<TARGET_FILE_DIR:unit_tests>
        COMMAND_EXPAND_LISTS
    )
endif()

# Add executable "test_dipole_operator"
add_executable(test_dipole_operator
    test_dipole_operator.cpp
)
target_link_libraries(test_dipole_operator PRIVATE pairinteraction)

if (WIN32)
    add_custom_command(TARGET test_dipole_operator POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:test_dipole_operator> $<TARGET_FILE_DIR:test_dipole_operator>
        COMMAND_EXPAND_LISTS
    )
endif()

# Add executable "test_system_atom"
add_executable(test_system_atom
    test_system_atom.cpp
)
target_link_libraries(test_system_atom PRIVATE pairinteraction)

if (WIN32)
    add_custom_command(TARGET test_system_atom POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:test_system_atom> $<TARGET_FILE_DIR:test_system_atom>
        COMMAND_EXPAND_LISTS
    )
endif()

# Add executable "test_starkmap"
add_executable(test_starkmap
test_starkmap.cpp
)
target_link_libraries(test_starkmap PRIVATE pairinteraction)

if (WIN32)
    add_custom_command(TARGET test_starkmap POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:test_starkmap> $<TARGET_FILE_DIR:test_starkmap>
        COMMAND_EXPAND_LISTS
    )
endif()

# Add tests building and calling the executables
include(ProcessorCount)
ProcessorCount(NPROC_DEFAULT)
set(NPROC ${NPROC_DEFAULT} CACHE STRING "Number of processors to use for building the tests")

add_test(NAME gen_unit_tests COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config "$<CONFIG>" --target unit_tests -j${NPROC})
add_test(NAME unit_tests  COMMAND unit_tests --database "${CMAKE_SOURCE_DIR}/data/database")
set_tests_properties(gen_unit_tests PROPERTIES FIXTURES_SETUP f_gen_unit_tests)
set_tests_properties(unit_tests PROPERTIES  FIXTURES_REQUIRED f_gen_unit_tests)

add_test(NAME gen_test_dipole_operator COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config "$<CONFIG>" --target test_dipole_operator -j${NPROC})
add_test(NAME test_dipole_operator  COMMAND test_dipole_operator --database "${CMAKE_SOURCE_DIR}/data/database")
set_tests_properties(gen_test_dipole_operator PROPERTIES FIXTURES_SETUP f_gen_test_dipole_operator)
set_tests_properties(test_dipole_operator PROPERTIES  FIXTURES_REQUIRED f_gen_test_dipole_operator)

add_test(NAME gen_test_system_atom COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config "$<CONFIG>" --target test_system_atom -j${NPROC})
add_test(NAME test_system_atom  COMMAND test_system_atom --database "${CMAKE_SOURCE_DIR}/data/database")
set_tests_properties(gen_test_system_atom PROPERTIES FIXTURES_SETUP f_gen_test_system_atom)
set_tests_properties(test_system_atom PROPERTIES  FIXTURES_REQUIRED f_gen_test_system_atom)

add_test(NAME gen_test_starkmap COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --config "$<CONFIG>" --target test_starkmap -j${NPROC})
add_test(NAME test_starkmap  COMMAND test_starkmap --database "${CMAKE_SOURCE_DIR}/data/database")
set_tests_properties(gen_test_starkmap PROPERTIES FIXTURES_SETUP f_gen_test_starkmap)
set_tests_properties(test_starkmap PROPERTIES  FIXTURES_REQUIRED f_gen_test_starkmap)

set_tests_properties(gen_unit_tests gen_test_dipole_operator gen_test_system_atom gen_test_starkmap PROPERTIES RESOURCE_LOCK compiler_access)
