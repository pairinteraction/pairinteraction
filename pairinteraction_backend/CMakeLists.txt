# Build the library
file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
list(FILTER SOURCES EXCLUDE REGEX ".*main\\.cpp$")

add_library(pairinteraction_backend SHARED ${SOURCES})

target_include_directories(pairinteraction_backend
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(pairinteraction_backend
    PRIVATE
        doctest::doctest
        Eigen3::Eigen
)

target_compile_definitions(pairinteraction_backend
    PRIVATE
        DOCTEST_CONFIG_SUPER_FAST_ASSERTS
        $<$<NOT:$<BOOL:${BUILD_TESTING}>>:DOCTEST_CONFIG_DISABLE>
)

target_compile_features(pairinteraction_backend
    PRIVATE
        cxx_std_20 # TODO is this fine or better -fconcepts?
)
set_target_properties(pairinteraction_backend
    PROPERTIES
        CXX_EXTENSIONS OFF
        CXX_CLANG_TIDY ""
        PREFIX ""
        INSTALL_RPATH $<$<BOOL:${SKBUILD}>:$<$<PLATFORM_ID:Linux>:$ORIGIN>$<$<PLATFORM_ID:Darwin>:@loader_path>>
)
target_include_directories(pairinteraction_backend
    PRIVATE
        "${Python_INCLUDE_DIRS};${Python_NumPy_INCLUDE_DIRS}"
)
target_compile_options(pairinteraction_backend
    PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/bigobj>
)
target_link_options(pairinteraction_backend
    PRIVATE
        $<$<CXX_COMPILER_ID:AppleClang>:-undefined dynamic_lookup>
        $<$<CXX_COMPILER_ID:Clang,GNU>:-Wl,--unresolved-symbols=ignore-all>
        $<$<CXX_COMPILER_ID:MSVC>:/NODEFAULTLIB:python3${Python_VERSION_MINOR}.lib /DEFAULTLIB:python3.lib /FORCE:UNRESOLVED>
)
target_link_libraries(pairinteraction_backend
    PRIVATE
        nanobind-static $<$<CXX_COMPILER_ID:MSVC>:${Python_LIBRARIES_GENERALIZED}>
)

if (BUILD_TESTING)

    # Build an executable for running the tests
    file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

    add_executable(pairinteraction_backend_test ${SOURCES})

    target_include_directories(pairinteraction_backend_test
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}
    )

    target_link_libraries(pairinteraction_backend_test
        PRIVATE
            doctest::doctest
            Eigen3::Eigen
    )
    target_compile_features(pairinteraction_backend_test
        PRIVATE
            cxx_std_20 # TODO is this fine or better -fconcepts?
    )
    set_target_properties(pairinteraction_backend_test
        PROPERTIES
            CXX_EXTENSIONS OFF
            CXX_CLANG_TIDY ""
            PREFIX ""
            INSTALL_RPATH $<$<BOOL:${SKBUILD}>:$<$<PLATFORM_ID:Linux>:$ORIGIN>$<$<PLATFORM_ID:Darwin>:@loader_path>>
    )
    target_include_directories(pairinteraction_backend_test
        PRIVATE
            "${Python_INCLUDE_DIRS};${Python_NumPy_INCLUDE_DIRS}"
    )
    target_compile_options(pairinteraction_backend_test
        PRIVATE
            $<$<CXX_COMPILER_ID:MSVC>:/bigobj>
    )
    target_link_options(pairinteraction_backend_test
        PRIVATE
            $<$<CXX_COMPILER_ID:AppleClang>:-undefined dynamic_lookup>
            $<$<CXX_COMPILER_ID:Clang,GNU>:-Wl,--unresolved-symbols=ignore-all>
            $<$<CXX_COMPILER_ID:MSVC>:/NODEFAULTLIB:python3${Python_VERSION_MINOR}.lib /DEFAULTLIB:python3.lib /FORCE:UNRESOLVED>
    )
    target_link_libraries(pairinteraction_backend_test
        PRIVATE
            nanobind-static $<$<CXX_COMPILER_ID:MSVC>:${Python_LIBRARIES_GENERALIZED}>
    )

    # TODO can we avoid building a separate executable but link the library to the test target?
    # add_executable(pairinteraction_backend_test main.cpp)
    # target_link_libraries(pairinteraction_backend_test
    #     PRIVATE
    #         pairinteraction_backend
    #         doctest::doctest
    # )

    # TODO can we use the following definitions for all targets that link to Eigen by changing the eigen target (also to include the MKL etc. target if necessary)
    # or should we set the definitions here?
    # EIGEN_USE_BLAS EIGEN_USE_LAPACKE
    # "-Dlapack_complex_float=std::complex<float>"
    # "-Dlapack_complex_double=std::complex<double>")
    # EIGEN_USE_MKL_ALL
    # "-DMKL_Complex8=std::complex<float>"
    # "-DMKL_Complex16=std::complex<double>")

    # Add tests
    # To only build and run these tests, use
    # "cmake --build . --target pairinteraction_backend_test && ctest --verbose -R pairinteraction_backend_test"
    add_test(NAME "pairinteraction_backend_test" COMMAND pairinteraction_backend_test)

endif()
