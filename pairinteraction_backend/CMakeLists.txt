# Add source files
file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
list(FILTER SOURCES EXCLUDE REGEX ".*main\\.cpp$")

# Build the library
add_library(pairinteraction_backend SHARED ${SOURCES})

target_include_directories(pairinteraction_backend
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(pairinteraction_backend
    PRIVATE
        doctest::doctest
)

target_compile_definitions(pairinteraction_backend
    PRIVATE
        DOCTEST_CONFIG_SUPER_FAST_ASSERTS
        $<$<NOT:$<BOOL:${BUILD_TESTING}>>:DOCTEST_CONFIG_DISABLE> 
)

if (BUILD_TESTING)

    # Build an executable for running the tests
    file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
    add_executable(pairinteraction_backend_test ${SOURCES})
    target_link_libraries(pairinteraction_backend_test
        PRIVATE
            doctest::doctest
    )

    # TODO can we avoid building a separate executable but link the library to the test target?
    # add_executable(pairinteraction_backend_test main.cpp)
    # target_link_libraries(pairinteraction_backend_test
    #     PRIVATE
    #         pairinteraction_backend
    #         doctest::doctest
    # )

    # Add tests
    # To only build and run these tests, use
    # "cmake --build . --target pairinteraction_backend_test && ctest --rerun-failed --output-on-failure -R pairinteraction_backend_test"
    add_test(NAME "pairinteraction_backend_test" COMMAND pairinteraction_backend_test)

endif()