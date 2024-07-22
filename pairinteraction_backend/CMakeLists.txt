set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Create backend object
add_library(backend OBJECT
  ./basis/Basis.cpp
  ./basis/Basis.hpp
  ./basis/BasisAtom.cpp
  ./basis/BasisAtom.hpp
  ./basis/BasisAtomCreator.cpp
  ./basis/BasisAtomCreator.hpp
  ./basis/BasisClassicalLight.cpp
  ./basis/BasisClassicalLight.hpp
  ./basis/BasisClassicalLightCreator.cpp
  ./basis/BasisClassicalLightCreator.hpp
  ./database/AtomDescriptionByParameters.hpp
  ./database/AtomDescriptionByRanges.hpp
  ./database/Database.cpp
  ./database/Database.hpp
  ./diagonalizer/Diagonalizer.cpp
  ./diagonalizer/Diagonalizer.hpp
  ./enums/OperatorType.hpp
  ./enums/Parity.cpp
  ./enums/Parity.hpp
  ./enums/TransformationType.hpp
  ./interfaces/TransformationBuilderInterface.cpp
  ./interfaces/TransformationBuilderInterface.hpp
  ./ket/Ket.cpp
  ./ket/Ket.hpp
  ./ket/KetAtom.cpp
  ./ket/KetAtom.hpp
  ./ket/KetAtomCreator.cpp
  ./ket/KetAtomCreator.hpp
  ./ket/KetClassicalLight.cpp
  ./ket/KetClassicalLight.hpp
  ./ket/KetClassicalLightCreator.cpp
  ./ket/KetClassicalLightCreator.hpp
  ./operator/Operator.cpp
  ./operator/Operator.hpp
  ./operator/OperatorAtom.cpp
  ./operator/OperatorAtom.hpp
  ./system/System.cpp
  ./system/System.hpp
  ./system/SystemAtom.cpp
  ./system/SystemAtom.hpp
  ./system/SystemClassicalLight.cpp
  ./system/SystemClassicalLight.hpp
  ./utils/args.hpp
  ./utils/eigen_assertion.hpp
  ./utils/eigen_compat.hpp
  ./utils/euler.cpp
  ./utils/euler.hpp
  ./utils/hash.hpp
  ./utils/ketid.hpp
  ./utils/maths.hpp
  ./utils/paths.hpp
  ./utils/Range.hpp
  ./utils/spherical.hpp
  ./utils/streamed.hpp
  ./utils/traits.hpp
  ./utils/wigner.cpp
  ./utils/wigner.hpp
  ./setup.cpp
  ./setup.hpp
  ./test.cpp
  ./test.hpp
)

target_link_libraries(backend PUBLIC
  $<$<BOOL:${WITH_COVERAGE}>:coverage::gcov>
  doctest::doctest
  duckdb::duckdb
  Eigen3::Eigen
  httplib::httplib
  $<$<TARGET_EXISTS:OpenSSL::SSL>:OpenSSL::SSL>
  $<$<TARGET_EXISTS:OpenSSL::Crypto>:OpenSSL::Crypto>
  $<$<TARGET_EXISTS:MKL::MKL>:MKL::MKL>
  nlohmann_json::nlohmann_json
  spdlog::spdlog
)

target_include_directories(backend PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_compile_definitions(backend PUBLIC
  DOCTEST_CONFIG_SUPER_FAST_ASSERTS
  $<$<TARGET_EXISTS:MKL::MKL>:EIGEN_USE_MKL_ALL>
  $<$<NOT:$<BOOL:${BUILD_TESTING}>>:DOCTEST_CONFIG_DISABLE>
)

target_compile_options(backend PUBLIC
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -Wpedantic>
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

target_compile_features(backend PUBLIC cxx_std_17)

# Create bindings object
add_library(bindings OBJECT
  # Headers
  ./basis/Basis.hpp
  ./basis/BasisAtom.hpp
  ./basis/BasisAtomCreator.hpp
  ./basis/BasisClassicalLight.hpp
  ./basis/BasisClassicalLightCreator.hpp
  ./database/AtomDescriptionByParameters.hpp
  ./database/AtomDescriptionByRanges.hpp
  ./database/Database.hpp
  ./diagonalizer/Diagonalizer.hpp
  ./enums/OperatorType.hpp
  ./enums/Parity.hpp
  ./enums/TransformationType.hpp
  ./interfaces/TransformationBuilderInterface.hpp
  ./ket/Ket.hpp
  ./ket/KetAtom.hpp
  ./ket/KetAtomCreator.hpp
  ./ket/KetClassicalLight.hpp
  ./ket/KetClassicalLightCreator.hpp
  ./operator/Operator.hpp
  ./operator/OperatorAtom.hpp
  ./system/System.hpp
  ./system/SystemAtom.hpp
  ./system/SystemClassicalLight.hpp
  ./utils/args.hpp
  ./utils/eigen_assertion.hpp
  ./utils/eigen_compat.hpp
  ./utils/euler.hpp
  ./utils/hash.hpp
  ./utils/ketid.hpp
  ./utils/maths.hpp
  ./utils/paths.hpp
  ./utils/Range.hpp
  ./utils/spherical.hpp
  ./utils/streamed.hpp
  ./utils/traits.hpp
  ./utils/wigner.hpp
  ./test.hpp
  ./setup.hpp
  # Bindings
  ./basis/Basis.py.cpp
  ./basis/Basis.py.hpp
  ./basis/BasisAtom.py.cpp
  ./basis/BasisAtom.py.hpp
  ./basis/BasisAtomCreator.py.cpp
  ./basis/BasisAtomCreator.py.hpp
  ./basis/BasisClassicalLight.py.cpp
  ./basis/BasisClassicalLight.py.hpp
  ./basis/BasisClassicalLightCreator.py.cpp
  ./basis/BasisClassicalLightCreator.py.hpp
  ./database/Database.py.cpp
  ./database/Database.py.hpp
  ./diagonalizer/Diagonalizer.py.cpp
  ./diagonalizer/Diagonalizer.py.hpp
  ./enums/OperatorType.py.cpp
  ./enums/OperatorType.py.hpp
  ./enums/TransformationType.py.cpp
  ./enums/TransformationType.py.hpp
  ./interfaces/TransformationBuilderInterface.py.cpp
  ./interfaces/TransformationBuilderInterface.py.hpp
  ./ket/Ket.py.cpp
  ./ket/Ket.py.hpp
  ./ket/KetAtom.py.cpp
  ./ket/KetAtom.py.hpp
  ./ket/KetAtomCreator.py.cpp
  ./ket/KetAtomCreator.py.hpp
  ./ket/KetClassicalLight.py.cpp
  ./ket/KetClassicalLight.py.hpp
  ./ket/KetClassicalLightCreator.py.cpp
  ./ket/KetClassicalLightCreator.py.hpp
  ./operator/Operator.py.cpp
  ./operator/Operator.py.hpp
  ./operator/OperatorAtom.py.cpp
  ./operator/OperatorAtom.py.hpp
  ./system/System.py.cpp
  ./system/System.py.hpp
  ./system/SystemAtom.py.cpp
  ./system/SystemAtom.py.hpp
  ./system/SystemClassicalLight.py.cpp
  ./system/SystemClassicalLight.py.hpp
  ./binding.cpp
  ./test.py.cpp
  ./test.py.hpp
)

target_include_directories(bindings PRIVATE "${Python_INCLUDE_DIRS};${Python_NumPy_INCLUDE_DIRS};${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(bindings PRIVATE $<$<BOOL:${WITH_COVERAGE}>:coverage::gcov> nanobind-static)
target_compile_options(bindings PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/bigobj /MP>)

# Assemble objects into libraries
add_library(pairinteraction_backend SHARED)

target_link_libraries(pairinteraction_backend
  PUBLIC $<$<CXX_COMPILER_ID:MSVC>:${Python_LIBRARIES_GENERALIZED}>
  PRIVATE backend bindings
)

set_target_properties(pairinteraction_backend PROPERTIES
  PREFIX ""
  SUFFIX "${NB_SUFFIX_EXT}"
  INSTALL_RPATH $<$<BOOL:${SKBUILD}>:$<$<PLATFORM_ID:Linux>:$ORIGIN>$<$<PLATFORM_ID:Darwin>:@loader_path>>
)

target_link_options(pairinteraction_backend PRIVATE
  $<$<CXX_COMPILER_ID:AppleClang>:-undefined dynamic_lookup>
  $<$<CXX_COMPILER_ID:Clang,GNU>:-Wl,--unresolved-symbols=ignore-all>
  $<$<CXX_COMPILER_ID:MSVC>:/NODEFAULTLIB:python3${Python_VERSION_MINOR}.lib /DEFAULTLIB:python3.lib /FORCE:UNRESOLVED>
)

if (WIN32)
  add_custom_command(TARGET pairinteraction_backend POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy -t $<TARGET_FILE_DIR:pairinteraction_backend> $<TARGET_RUNTIME_DLLS:pairinteraction_backend>
    COMMAND_EXPAND_LISTS
  )
endif()

install(IMPORTED_RUNTIME_ARTIFACTS duckdb DESTINATION pairinteraction_backend)
install(TARGETS pairinteraction_backend LIBRARY DESTINATION pairinteraction_backend)

# Build an executable for running the tests
if (BUILD_TESTING)
  add_executable(pairinteraction_backend_test
    main.cpp
  )
  target_link_libraries(pairinteraction_backend_test PRIVATE backend)

  if (WIN32)
    add_custom_command(TARGET pairinteraction_backend_test POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy -t $<TARGET_FILE_DIR:pairinteraction_backend_test> $<TARGET_RUNTIME_DLLS:pairinteraction_backend_test>
      COMMAND_EXPAND_LISTS
    )
  endif()

  # Add tests
  add_test(NAME "pairinteraction_backend_test" COMMAND pairinteraction_backend_test --database "${CMAKE_SOURCE_DIR}/test/data/database")
endif()
