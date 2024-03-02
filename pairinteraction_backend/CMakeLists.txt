# Build object files for the backend
file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
list(FILTER SOURCES EXCLUDE REGEX ".*main\\.cpp$")

add_library(pairinteraction_backend OBJECT ${SOURCES})

target_include_directories(pairinteraction_backend
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        "${Python_INCLUDE_DIRS};${Python_NumPy_INCLUDE_DIRS}"
)

target_link_libraries(pairinteraction_backend
    PUBLIC
        doctest::doctest
        Eigen3::Eigen
        $<$<CXX_COMPILER_ID:MSVC>:${Python_LIBRARIES_GENERALIZED}>
    PRIVATE
        nanobind-static
)

target_compile_definitions(pairinteraction_backend
    PUBLIC
        DOCTEST_CONFIG_SUPER_FAST_ASSERTS
        $<$<NOT:$<BOOL:${BUILD_TESTING}>>:DOCTEST_CONFIG_DISABLE>
)

target_compile_features(pairinteraction_backend
    PUBLIC
        cxx_std_17
)

set_target_properties(pairinteraction_backend
    PROPERTIES
        CXX_EXTENSIONS OFF
        CXX_CLANG_TIDY ""
        PREFIX ""
        INSTALL_RPATH $<$<BOOL:${SKBUILD}>:$<$<PLATFORM_ID:Linux>:$ORIGIN>$<$<PLATFORM_ID:Darwin>:@loader_path>>
)

target_compile_options(pairinteraction_backend
    PUBLIC
        $<$<CXX_COMPILER_ID:MSVC>:/bigobj>
        -fPIC
)

target_link_options(pairinteraction_backend
    PUBLIC
        $<$<CXX_COMPILER_ID:AppleClang>:-undefined dynamic_lookup>
        $<$<CXX_COMPILER_ID:Clang,GNU>:-Wl,--unresolved-symbols=ignore-all>
        $<$<CXX_COMPILER_ID:MSVC>:/NODEFAULTLIB:python3${Python_VERSION_MINOR}.lib /DEFAULTLIB:python3.lib /FORCE:UNRESOLVED>
)

# Build a shared library
add_library(pairinteraction_backend_shared SHARED $<TARGET_OBJECTS:pairinteraction_backend>)
target_link_libraries(pairinteraction_backend_shared
PRIVATE
    nanobind-static
)

# Build an executable for running the tests
if (BUILD_TESTING)
    add_executable(pairinteraction_backend_test main.cpp
        $<TARGET_OBJECTS:pairinteraction_backend>
    )

    target_link_libraries(pairinteraction_backend_test
        PRIVATE
            nanobind-static
    )

    # Add tests
    # To only build and run these tests, use
    # "cmake --build . --target pairinteraction_backend_test && ctest --verbose -R pairinteraction_backend_test"
    add_test(NAME "pairinteraction_backend_test" COMMAND pairinteraction_backend_test)

endif()
