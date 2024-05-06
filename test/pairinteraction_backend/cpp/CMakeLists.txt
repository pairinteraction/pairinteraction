# Build executable "test_dipole_operator"
add_executable(test_dipole_operator
    $<TARGET_OBJECTS:backend>
    test_dipole_operator.cpp
)

target_compile_definitions(test_dipole_operator PRIVATE
    DOCTEST_CONFIG_SUPER_FAST_ASSERTS
    JSON_USE_IMPLICIT_CONVERSIONS=0
)

target_link_libraries(test_dipole_operator
    PRIVATE spdlog::spdlog duckdb::duckdb httplib::httplib Eigen3::Eigen
)

target_compile_options(test_dipole_operator PRIVATE
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

target_include_directories(test_dipole_operator PRIVATE
  "${CMAKE_SOURCE_DIR}/pairinteraction_backend"
)

if (WIN32)
    add_custom_command(TARGET test_dipole_operator POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy -t $<TARGET_FILE_DIR:test_dipole_operator> $<TARGET_RUNTIME_DLLS:test_dipole_operator>
        COMMAND_EXPAND_LISTS
    )
endif()

target_compile_features(test_dipole_operator PUBLIC cxx_std_17)

# Add tests calling the executables
add_test(NAME "test_dipole_operator"  COMMAND test_dipole_operator)
