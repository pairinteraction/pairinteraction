# Build executable "test_dipole_operator"
add_executable(test_dipole_operator
    test_dipole_operator.cpp
)
target_link_libraries(test_dipole_operator PRIVATE backend)
target_compile_definitions(test_dipole_operator PRIVATE
    PAIRINTERACTION_LOCAL_DATABASE_DIR="${CMAKE_SOURCE_DIR}/test/database"
)

if (WIN32)
    add_custom_command(TARGET test_dipole_operator POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy -t $<TARGET_FILE_DIR:test_dipole_operator> $<TARGET_RUNTIME_DLLS:test_dipole_operator>
        COMMAND_EXPAND_LISTS
    )
endif()

# Add tests calling the executables
add_test(NAME "test_dipole_operator"  COMMAND test_dipole_operator)
