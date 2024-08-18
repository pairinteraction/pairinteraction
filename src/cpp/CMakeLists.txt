set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Create pairinteraction library object
add_library(pairinteraction OBJECT
  ./pairinteraction/basis/Basis.cpp
  ./pairinteraction/basis/Basis.hpp
  ./pairinteraction/basis/BasisAtom.cpp
  ./pairinteraction/basis/BasisAtom.hpp
  ./pairinteraction/basis/BasisAtomCreator.cpp
  ./pairinteraction/basis/BasisAtomCreator.hpp
  ./pairinteraction/basis/BasisAtomCreator.test.cpp
  ./pairinteraction/basis/BasisClassicalLight.cpp
  ./pairinteraction/basis/BasisClassicalLight.hpp
  ./pairinteraction/basis/BasisClassicalLightCreator.cpp
  ./pairinteraction/basis/BasisClassicalLightCreator.hpp
  ./pairinteraction/basis/BasisClassicalLightCreator.test.cpp
  ./pairinteraction/database/AtomDescriptionByParameters.hpp
  ./pairinteraction/database/AtomDescriptionByRanges.hpp
  ./pairinteraction/database/Database.cpp
  ./pairinteraction/database/Database.hpp
  ./pairinteraction/database/Database.test.cpp
  ./pairinteraction/diagonalizer/diagonalize.cpp
  ./pairinteraction/diagonalizer/diagonalize.hpp
  ./pairinteraction/diagonalizer/DiagonalizerEigen.cpp
  ./pairinteraction/diagonalizer/DiagonalizerEigen.hpp
  ./pairinteraction/diagonalizer/DiagonalizerFeast.cpp
  ./pairinteraction/diagonalizer/DiagonalizerFeast.hpp
  ./pairinteraction/diagonalizer/DiagonalizerLapacke.cpp
  ./pairinteraction/diagonalizer/DiagonalizerLapacke.hpp
  ./pairinteraction/enums/OperatorType.hpp
  ./pairinteraction/enums/Parity.hpp
  ./pairinteraction/enums/Parity.test.cpp
  ./pairinteraction/enums/TransformationType.hpp
  ./pairinteraction/interfaces/DiagonalizerInterface.cpp
  ./pairinteraction/interfaces/DiagonalizerInterface.hpp
  ./pairinteraction/interfaces/TransformationBuilderInterface.cpp
  ./pairinteraction/interfaces/TransformationBuilderInterface.hpp
  ./pairinteraction/ket/Ket.cpp
  ./pairinteraction/ket/Ket.hpp
  ./pairinteraction/ket/Ket.test.cpp
  ./pairinteraction/ket/KetAtom.cpp
  ./pairinteraction/ket/KetAtom.hpp
  ./pairinteraction/ket/KetAtomCreator.cpp
  ./pairinteraction/ket/KetAtomCreator.hpp
  ./pairinteraction/ket/KetAtomCreator.test.cpp
  ./pairinteraction/ket/KetClassicalLight.cpp
  ./pairinteraction/ket/KetClassicalLight.hpp
  ./pairinteraction/ket/KetClassicalLightCreator.cpp
  ./pairinteraction/ket/KetClassicalLightCreator.hpp
  ./pairinteraction/ket/KetClassicalLightCreator.test.cpp
  ./pairinteraction/operator/Operator.cpp
  ./pairinteraction/operator/Operator.hpp
  ./pairinteraction/operator/OperatorAtom.cpp
  ./pairinteraction/operator/OperatorAtom.hpp
  ./pairinteraction/system/System.cpp
  ./pairinteraction/system/System.hpp
  ./pairinteraction/system/SystemAtom.cpp
  ./pairinteraction/system/SystemAtom.hpp
  ./pairinteraction/system/SystemAtom.test.cpp
  ./pairinteraction/system/SystemClassicalLight.cpp
  ./pairinteraction/system/SystemClassicalLight.hpp
  ./pairinteraction/tools/setup.cpp
  ./pairinteraction/tools/setup.hpp
  ./pairinteraction/tools/test.cpp
  ./pairinteraction/tools/test.hpp
  ./pairinteraction/utils/args.hpp
  ./pairinteraction/utils/eigen_assertion.hpp
  ./pairinteraction/utils/eigen_compat.hpp
  ./pairinteraction/utils/euler.hpp
  ./pairinteraction/utils/euler.test.cpp
  ./pairinteraction/utils/hash.hpp
  ./pairinteraction/utils/ketid.hpp
  ./pairinteraction/utils/maths.hpp
  ./pairinteraction/utils/paths.hpp
  ./pairinteraction/utils/Range.hpp
  ./pairinteraction/utils/spherical.hpp
  ./pairinteraction/utils/streamed.hpp
  ./pairinteraction/utils/traits.hpp
  ./pairinteraction/utils/wigner.hpp
  ./pairinteraction/utils/wigner.test.cpp
  ./pairinteraction/pairinteraction.hpp
)

target_link_libraries(pairinteraction PUBLIC
  $<$<BOOL:${WITH_COVERAGE}>:coverage::gcov>
  doctest::doctest
  duckdb::duckdb
  Eigen3::Eigen
  httplib::httplib
  nlohmann_json::nlohmann_json
  spdlog::spdlog
  TBB::tbb
  $<$<TARGET_EXISTS:MKL::MKL>:MKL::MKL>
  $<$<AND:$<TARGET_EXISTS:LAPACKE::LAPACKE>,$<NOT:$<TARGET_EXISTS:MKL::MKL>>>:LAPACKE::LAPACKE>
  $<$<TARGET_EXISTS:OpenSSL::SSL>:OpenSSL::SSL>
  $<$<TARGET_EXISTS:OpenSSL::Crypto>:OpenSSL::Crypto>
)

target_include_directories(pairinteraction PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_compile_definitions(pairinteraction PUBLIC
  DOCTEST_CONFIG_SUPER_FAST_ASSERTS
  $<$<TARGET_EXISTS:MKL::MKL>:WITH_MKL WITH_LAPACKE EIGEN_USE_MKL_ALL EIGEN_MKL_NO_DIRECT_CALL>
  $<$<AND:$<TARGET_EXISTS:LAPACKE::LAPACKE>,$<NOT:$<TARGET_EXISTS:MKL::MKL>>>:WITH_LAPACKE EIGEN_USE_BLAS EIGEN_USE_LAPACKE>
  $<$<NOT:$<BOOL:${BUILD_TESTING}>>:DOCTEST_CONFIG_DISABLE>
)

if (TARGET MKL::MKL)
  target_compile_definitions(pairinteraction PUBLIC
    MKL_Complex8=std::complex<float>
    MKL_Complex16=std::complex<double>
  )
elseif (TARGET LAPACKE::LAPACKE)
  target_compile_definitions(pairinteraction PUBLIC
    lapack_complex_float=std::complex<float>
    lapack_complex_double=std::complex<double>
  )
endif()

target_compile_options(pairinteraction PUBLIC
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -Wpedantic>
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

target_compile_features(pairinteraction PUBLIC cxx_std_17)

# Create bindings object
add_library(bindings OBJECT
  # Headers
  ./pairinteraction/basis/Basis.hpp
  ./pairinteraction/basis/BasisAtom.hpp
  ./pairinteraction/basis/BasisAtomCreator.hpp
  ./pairinteraction/basis/BasisClassicalLight.hpp
  ./pairinteraction/basis/BasisClassicalLightCreator.hpp
  ./pairinteraction/database/AtomDescriptionByParameters.hpp
  ./pairinteraction/database/AtomDescriptionByRanges.hpp
  ./pairinteraction/database/Database.hpp
  ./pairinteraction/diagonalizer/diagonalize.hpp
  ./pairinteraction/diagonalizer/DiagonalizerEigen.hpp
  ./pairinteraction/diagonalizer/DiagonalizerFeast.hpp
  ./pairinteraction/diagonalizer/DiagonalizerLapacke.hpp
  ./pairinteraction/enums/OperatorType.hpp
  ./pairinteraction/enums/Parity.hpp
  ./pairinteraction/enums/TransformationType.hpp
  ./pairinteraction/interfaces/DiagonalizerInterface.hpp
  ./pairinteraction/interfaces/TransformationBuilderInterface.hpp
  ./pairinteraction/ket/Ket.hpp
  ./pairinteraction/ket/KetAtom.hpp
  ./pairinteraction/ket/KetAtomCreator.hpp
  ./pairinteraction/ket/KetClassicalLight.hpp
  ./pairinteraction/ket/KetClassicalLightCreator.hpp
  ./pairinteraction/operator/Operator.hpp
  ./pairinteraction/operator/OperatorAtom.hpp
  ./pairinteraction/system/System.hpp
  ./pairinteraction/system/SystemAtom.hpp
  ./pairinteraction/system/SystemClassicalLight.hpp
  ./pairinteraction/tools/test.hpp
  ./pairinteraction/tools/setup.hpp
  ./pairinteraction/utils/args.hpp
  ./pairinteraction/utils/eigen_assertion.hpp
  ./pairinteraction/utils/eigen_compat.hpp
  ./pairinteraction/utils/euler.hpp
  ./pairinteraction/utils/hash.hpp
  ./pairinteraction/utils/ketid.hpp
  ./pairinteraction/utils/maths.hpp
  ./pairinteraction/utils/paths.hpp
  ./pairinteraction/utils/Range.hpp
  ./pairinteraction/utils/spherical.hpp
  ./pairinteraction/utils/streamed.hpp
  ./pairinteraction/utils/traits.hpp
  ./pairinteraction/utils/wigner.hpp
  # Bindings
  ./pairinteraction/basis/Basis.py.cpp
  ./pairinteraction/basis/Basis.py.hpp
  ./pairinteraction/basis/BasisAtom.py.cpp
  ./pairinteraction/basis/BasisAtom.py.hpp
  ./pairinteraction/basis/BasisAtomCreator.py.cpp
  ./pairinteraction/basis/BasisAtomCreator.py.hpp
  ./pairinteraction/basis/BasisClassicalLight.py.cpp
  ./pairinteraction/basis/BasisClassicalLight.py.hpp
  ./pairinteraction/basis/BasisClassicalLightCreator.py.cpp
  ./pairinteraction/basis/BasisClassicalLightCreator.py.hpp
  ./pairinteraction/database/Database.py.cpp
  ./pairinteraction/database/Database.py.hpp
  ./pairinteraction/enums/OperatorType.py.cpp
  ./pairinteraction/enums/OperatorType.py.hpp
  ./pairinteraction/enums/TransformationType.py.cpp
  ./pairinteraction/enums/TransformationType.py.hpp
  ./pairinteraction/interfaces/TransformationBuilderInterface.py.cpp
  ./pairinteraction/interfaces/TransformationBuilderInterface.py.hpp
  ./pairinteraction/ket/Ket.py.cpp
  ./pairinteraction/ket/Ket.py.hpp
  ./pairinteraction/ket/KetAtom.py.cpp
  ./pairinteraction/ket/KetAtom.py.hpp
  ./pairinteraction/ket/KetAtomCreator.py.cpp
  ./pairinteraction/ket/KetAtomCreator.py.hpp
  ./pairinteraction/ket/KetClassicalLight.py.cpp
  ./pairinteraction/ket/KetClassicalLight.py.hpp
  ./pairinteraction/ket/KetClassicalLightCreator.py.cpp
  ./pairinteraction/ket/KetClassicalLightCreator.py.hpp
  ./pairinteraction/operator/Operator.py.cpp
  ./pairinteraction/operator/Operator.py.hpp
  ./pairinteraction/operator/OperatorAtom.py.cpp
  ./pairinteraction/operator/OperatorAtom.py.hpp
  ./pairinteraction/system/System.py.cpp
  ./pairinteraction/system/System.py.hpp
  ./pairinteraction/system/SystemAtom.py.cpp
  ./pairinteraction/system/SystemAtom.py.hpp
  ./pairinteraction/system/SystemClassicalLight.py.cpp
  ./pairinteraction/system/SystemClassicalLight.py.hpp
  ./pairinteraction/tools/test.py.cpp
  ./pairinteraction/tools/test.py.hpp
  ./pairinteraction/binding.cpp
)

target_include_directories(bindings PRIVATE "${Python_INCLUDE_DIRS};${Python_NumPy_INCLUDE_DIRS};${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(bindings PRIVATE $<$<BOOL:${WITH_COVERAGE}>:coverage::gcov> nanobind-static)
target_compile_options(bindings PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/bigobj /MP>)

# Assemble objects into backend library
add_library(backend SHARED)

target_link_libraries(backend
  PUBLIC $<$<CXX_COMPILER_ID:MSVC>:${Python_LIBRARIES_GENERALIZED}>
  PRIVATE pairinteraction bindings
)

set_target_properties(backend PROPERTIES
  PREFIX ""
  SUFFIX "${NB_SUFFIX_EXT}"
  INSTALL_RPATH $<$<BOOL:${SKBUILD}>:$<$<PLATFORM_ID:Linux>:$ORIGIN>$<$<PLATFORM_ID:Darwin>:@loader_path>>
)

target_link_options(backend PRIVATE
  $<$<CXX_COMPILER_ID:AppleClang>:-undefined dynamic_lookup>
  $<$<CXX_COMPILER_ID:Clang,GNU>:-Wl,--unresolved-symbols=ignore-all>
  $<$<CXX_COMPILER_ID:MSVC>:/NODEFAULTLIB:python3${Python_VERSION_MINOR}.lib /DEFAULTLIB:python3.lib /FORCE:UNRESOLVED>
)

if (WIN32)
  add_custom_command(TARGET backend POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy -t $<TARGET_FILE_DIR:backend> $<TARGET_RUNTIME_DLLS:backend>
    COMMAND_EXPAND_LISTS
  )
endif()

install(IMPORTED_RUNTIME_ARTIFACTS duckdb DESTINATION pairinteraction)
install(TARGETS backend DESTINATION pairinteraction)
if (WIN32)
  install(FILES $<TARGET_RUNTIME_DLLS:backend> DESTINATION pairinteraction)
endif()

# Build an executable for running the tests
if (BUILD_TESTING)
  add_executable(pairinteraction_test
    ./pairinteraction_test/main.cpp
  )
  target_link_libraries(pairinteraction_test PRIVATE pairinteraction)

  if (WIN32)
    add_custom_command(TARGET pairinteraction_test POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy -t $<TARGET_FILE_DIR:pairinteraction_test> $<TARGET_RUNTIME_DLLS:pairinteraction_test>
      COMMAND_EXPAND_LISTS
    )
  endif()

  # Add tests
  add_test(NAME "pairinteraction_test" COMMAND pairinteraction_test --database "${CMAKE_SOURCE_DIR}/test/data/database")
endif()
