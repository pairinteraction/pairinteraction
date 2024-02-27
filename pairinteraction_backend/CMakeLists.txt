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
