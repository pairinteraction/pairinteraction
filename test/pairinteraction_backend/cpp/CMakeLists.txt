# Build executable "test_dipole_operator"
add_executable(test_dipole_operator
    test_dipole_operator.cpp
)
target_link_libraries(test_dipole_operator PRIVATE backend)

if (WIN32)
    add_custom_command(TARGET test_dipole_operator POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy -t $<TARGET_FILE_DIR:test_dipole_operator> $<TARGET_RUNTIME_DLLS:test_dipole_operator>
        COMMAND_EXPAND_LISTS
    )
endif()

# Build executable "system_atom"
add_executable(test_system_atom
    test_system_atom.cpp
)
target_link_libraries(test_system_atom PRIVATE backend)

if (WIN32)
    add_custom_command(TARGET test_system_atom POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy -t $<TARGET_FILE_DIR:test_system_atom> $<TARGET_RUNTIME_DLLS:test_system_atom>
        COMMAND_EXPAND_LISTS
    )
endif()

# Add tests calling the executables
add_test(NAME "test_dipole_operator"  COMMAND test_dipole_operator --database "${CMAKE_SOURCE_DIR}/test/data/database")
add_test(NAME "test_system_atom"  COMMAND test_system_atom --database "${CMAKE_SOURCE_DIR}/test/data/database")
