set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Create pairinteraction library object
add_library(pintr OBJECT
  ./pintr/basis/Basis.cpp
  ./pintr/basis/Basis.hpp
  ./pintr/basis/BasisAtom.cpp
  ./pintr/basis/BasisAtom.hpp
  ./pintr/basis/BasisAtomCreator.cpp
  ./pintr/basis/BasisAtomCreator.hpp
  ./pintr/basis/BasisAtomCreator.test.cpp
  ./pintr/basis/BasisClassicalLight.cpp
  ./pintr/basis/BasisClassicalLight.hpp
  ./pintr/basis/BasisClassicalLightCreator.cpp
  ./pintr/basis/BasisClassicalLightCreator.hpp
  ./pintr/basis/BasisClassicalLightCreator.test.cpp
  ./pintr/database/AtomDescriptionByParameters.hpp
  ./pintr/database/AtomDescriptionByRanges.hpp
  ./pintr/database/Database.cpp
  ./pintr/database/Database.hpp
  ./pintr/database/Database.test.cpp
  ./pintr/diagonalizer/diagonalize.cpp
  ./pintr/diagonalizer/diagonalize.hpp
  ./pintr/diagonalizer/DiagonalizerEigen.cpp
  ./pintr/diagonalizer/DiagonalizerEigen.hpp
  ./pintr/diagonalizer/DiagonalizerFeast.cpp
  ./pintr/diagonalizer/DiagonalizerFeast.hpp
  ./pintr/diagonalizer/DiagonalizerLapacke.cpp
  ./pintr/diagonalizer/DiagonalizerLapacke.hpp
  ./pintr/enums/OperatorType.hpp
  ./pintr/enums/Parity.hpp
  ./pintr/enums/Parity.test.cpp
  ./pintr/enums/TransformationType.hpp
  ./pintr/interfaces/DiagonalizerInterface.cpp
  ./pintr/interfaces/DiagonalizerInterface.hpp
  ./pintr/interfaces/TransformationBuilderInterface.cpp
  ./pintr/interfaces/TransformationBuilderInterface.hpp
  ./pintr/ket/Ket.cpp
  ./pintr/ket/Ket.hpp
  ./pintr/ket/Ket.test.cpp
  ./pintr/ket/KetAtom.cpp
  ./pintr/ket/KetAtom.hpp
  ./pintr/ket/KetAtomCreator.cpp
  ./pintr/ket/KetAtomCreator.hpp
  ./pintr/ket/KetAtomCreator.test.cpp
  ./pintr/ket/KetClassicalLight.cpp
  ./pintr/ket/KetClassicalLight.hpp
  ./pintr/ket/KetClassicalLightCreator.cpp
  ./pintr/ket/KetClassicalLightCreator.hpp
  ./pintr/ket/KetClassicalLightCreator.test.cpp
  ./pintr/operator/Operator.cpp
  ./pintr/operator/Operator.hpp
  ./pintr/operator/OperatorAtom.cpp
  ./pintr/operator/OperatorAtom.hpp
  ./pintr/system/System.cpp
  ./pintr/system/System.hpp
  ./pintr/system/SystemAtom.cpp
  ./pintr/system/SystemAtom.hpp
  ./pintr/system/SystemAtom.test.cpp
  ./pintr/system/SystemClassicalLight.cpp
  ./pintr/system/SystemClassicalLight.hpp
  ./pintr/tools/setup.cpp
  ./pintr/tools/setup.hpp
  ./pintr/tools/test.cpp
  ./pintr/tools/test.hpp
  ./pintr/utils/args.hpp
  ./pintr/utils/eigen_assertion.hpp
  ./pintr/utils/eigen_compat.hpp
  ./pintr/utils/euler.hpp
  ./pintr/utils/euler.test.cpp
  ./pintr/utils/hash.hpp
  ./pintr/utils/ketid.hpp
  ./pintr/utils/maths.hpp
  ./pintr/utils/paths.hpp
  ./pintr/utils/Range.hpp
  ./pintr/utils/spherical.hpp
  ./pintr/utils/streamed.hpp
  ./pintr/utils/traits.hpp
  ./pintr/utils/wigner.hpp
  ./pintr/utils/wigner.test.cpp
  ./pintr/pintr.hpp
)

target_link_libraries(pintr PUBLIC
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

target_include_directories(pintr PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_compile_definitions(pintr PUBLIC
  DOCTEST_CONFIG_SUPER_FAST_ASSERTS
  $<$<TARGET_EXISTS:MKL::MKL>:WITH_MKL WITH_LAPACKE EIGEN_USE_MKL_ALL EIGEN_MKL_NO_DIRECT_CALL>
  $<$<AND:$<TARGET_EXISTS:LAPACKE::LAPACKE>,$<NOT:$<TARGET_EXISTS:MKL::MKL>>>:WITH_LAPACKE EIGEN_USE_BLAS EIGEN_USE_LAPACKE>
  $<$<NOT:$<BOOL:${BUILD_TESTING}>>:DOCTEST_CONFIG_DISABLE>
)

if (TARGET MKL::MKL)
  target_compile_definitions(pintr PUBLIC
    MKL_Complex8=std::complex<float>
    MKL_Complex16=std::complex<double>
  )
elseif (TARGET LAPACKE::LAPACKE)
  target_compile_definitions(pintr PUBLIC
    lapack_complex_float=std::complex<float>
    lapack_complex_double=std::complex<double>
  )
endif()

target_compile_options(pintr PUBLIC
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -Wpedantic>
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

target_compile_features(pintr PUBLIC cxx_std_17)

# Create bindings object
add_library(bindings OBJECT
  # Headers
  ./pintr/basis/Basis.hpp
  ./pintr/basis/BasisAtom.hpp
  ./pintr/basis/BasisAtomCreator.hpp
  ./pintr/basis/BasisClassicalLight.hpp
  ./pintr/basis/BasisClassicalLightCreator.hpp
  ./pintr/database/AtomDescriptionByParameters.hpp
  ./pintr/database/AtomDescriptionByRanges.hpp
  ./pintr/database/Database.hpp
  ./pintr/diagonalizer/diagonalize.hpp
  ./pintr/diagonalizer/DiagonalizerEigen.hpp
  ./pintr/diagonalizer/DiagonalizerFeast.hpp
  ./pintr/diagonalizer/DiagonalizerLapacke.hpp
  ./pintr/enums/OperatorType.hpp
  ./pintr/enums/Parity.hpp
  ./pintr/enums/TransformationType.hpp
  ./pintr/interfaces/DiagonalizerInterface.hpp
  ./pintr/interfaces/TransformationBuilderInterface.hpp
  ./pintr/ket/Ket.hpp
  ./pintr/ket/KetAtom.hpp
  ./pintr/ket/KetAtomCreator.hpp
  ./pintr/ket/KetClassicalLight.hpp
  ./pintr/ket/KetClassicalLightCreator.hpp
  ./pintr/operator/Operator.hpp
  ./pintr/operator/OperatorAtom.hpp
  ./pintr/system/System.hpp
  ./pintr/system/SystemAtom.hpp
  ./pintr/system/SystemClassicalLight.hpp
  ./pintr/tools/test.hpp
  ./pintr/tools/setup.hpp
  ./pintr/utils/args.hpp
  ./pintr/utils/eigen_assertion.hpp
  ./pintr/utils/eigen_compat.hpp
  ./pintr/utils/euler.hpp
  ./pintr/utils/hash.hpp
  ./pintr/utils/ketid.hpp
  ./pintr/utils/maths.hpp
  ./pintr/utils/paths.hpp
  ./pintr/utils/Range.hpp
  ./pintr/utils/spherical.hpp
  ./pintr/utils/streamed.hpp
  ./pintr/utils/traits.hpp
  ./pintr/utils/wigner.hpp
  # Bindings
  ./pintr/basis/Basis.py.cpp
  ./pintr/basis/Basis.py.hpp
  ./pintr/basis/BasisAtom.py.cpp
  ./pintr/basis/BasisAtom.py.hpp
  ./pintr/basis/BasisAtomCreator.py.cpp
  ./pintr/basis/BasisAtomCreator.py.hpp
  ./pintr/basis/BasisClassicalLight.py.cpp
  ./pintr/basis/BasisClassicalLight.py.hpp
  ./pintr/basis/BasisClassicalLightCreator.py.cpp
  ./pintr/basis/BasisClassicalLightCreator.py.hpp
  ./pintr/database/Database.py.cpp
  ./pintr/database/Database.py.hpp
  ./pintr/enums/OperatorType.py.cpp
  ./pintr/enums/OperatorType.py.hpp
  ./pintr/enums/TransformationType.py.cpp
  ./pintr/enums/TransformationType.py.hpp
  ./pintr/interfaces/TransformationBuilderInterface.py.cpp
  ./pintr/interfaces/TransformationBuilderInterface.py.hpp
  ./pintr/ket/Ket.py.cpp
  ./pintr/ket/Ket.py.hpp
  ./pintr/ket/KetAtom.py.cpp
  ./pintr/ket/KetAtom.py.hpp
  ./pintr/ket/KetAtomCreator.py.cpp
  ./pintr/ket/KetAtomCreator.py.hpp
  ./pintr/ket/KetClassicalLight.py.cpp
  ./pintr/ket/KetClassicalLight.py.hpp
  ./pintr/ket/KetClassicalLightCreator.py.cpp
  ./pintr/ket/KetClassicalLightCreator.py.hpp
  ./pintr/operator/Operator.py.cpp
  ./pintr/operator/Operator.py.hpp
  ./pintr/operator/OperatorAtom.py.cpp
  ./pintr/operator/OperatorAtom.py.hpp
  ./pintr/system/System.py.cpp
  ./pintr/system/System.py.hpp
  ./pintr/system/SystemAtom.py.cpp
  ./pintr/system/SystemAtom.py.hpp
  ./pintr/system/SystemClassicalLight.py.cpp
  ./pintr/system/SystemClassicalLight.py.hpp
  ./pintr/tools/test.py.cpp
  ./pintr/tools/test.py.hpp
  ./pintr/binding.cpp
)

target_include_directories(bindings PRIVATE "${Python_INCLUDE_DIRS};${Python_NumPy_INCLUDE_DIRS};${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(bindings PRIVATE $<$<BOOL:${WITH_COVERAGE}>:coverage::gcov> nanobind-static)
target_compile_options(bindings PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/bigobj /MP>)

# Assemble objects into backend library
add_library(backend SHARED)

target_link_libraries(backend
  PUBLIC $<$<CXX_COMPILER_ID:MSVC>:${Python_LIBRARIES_GENERALIZED}>
  PRIVATE pintr bindings
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
  add_executable(pintr_test
    ./pintr_test/main.cpp
  )
  target_link_libraries(pintr_test PRIVATE pintr)

  if (WIN32)
    add_custom_command(TARGET pintr_test POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy -t $<TARGET_FILE_DIR:pintr_test> $<TARGET_RUNTIME_DLLS:pintr_test>
      COMMAND_EXPAND_LISTS
    )
  endif()

  # Add tests
  add_test(NAME "pintr_test" COMMAND pintr_test --database "${CMAKE_SOURCE_DIR}/test/data/database")
endif()
