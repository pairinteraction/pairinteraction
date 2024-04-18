# Create backend object
add_library(backend OBJECT)
target_sources(backend PUBLIC
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
  ./cache/Cache.cpp
  ./cache/Cache.hpp
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
  ./system/System.cpp
  ./system/System.hpp
  ./system/SystemAtom.cpp
  ./system/SystemAtom.hpp
  ./system/SystemClassicalLight.cpp
  ./system/SystemClassicalLight.hpp
  ./utils/euler.hpp
  ./utils/maths.hpp
  ./utils/traits.hpp
  ./utils/wigner.hpp
)

target_include_directories(backend PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(backend
  PUBLIC doctest::doctest Eigen3::Eigen
  PRIVATE spdlog::spdlog
)

target_compile_definitions(backend PUBLIC
  DOCTEST_CONFIG_SUPER_FAST_ASSERTS
  $<$<NOT:$<BOOL:${BUILD_TESTING}>>:DOCTEST_CONFIG_DISABLE>
)

target_compile_features(backend PUBLIC cxx_std_20)

set_target_properties(backend PROPERTIES
  CXX_EXTENSIONS OFF
  POSITION_INDEPENDENT_CODE ON
)

# Create bindings object
add_library(bindings OBJECT)
target_sources(bindings PUBLIC
  # Headers
  ./basis/Basis.hpp
  ./basis/BasisAtom.hpp
  ./basis/BasisAtomCreator.hpp
  ./basis/BasisClassicalLight.hpp
  ./basis/BasisClassicalLightCreator.hpp
  ./cache/Cache.cpp
  ./cache/Cache.hpp
  ./ket/Ket.hpp
  ./ket/KetAtom.hpp
  ./ket/KetAtomCreator.hpp
  ./ket/KetClassicalLight.hpp
  ./ket/KetClassicalLightCreator.hpp
  ./system/System.hpp
  ./system/SystemAtom.hpp
  ./system/SystemClassicalLight.hpp
  ./utils/euler.hpp
  ./utils/maths.hpp
  ./utils/traits.hpp
  ./utils/wigner.hpp
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
  ./cache/Cache.py.cpp
  ./cache/Cache.py.hpp
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
  ./system/System.py.cpp
  ./system/System.py.hpp
  ./system/SystemAtom.py.cpp
  ./system/SystemAtom.py.hpp
  ./system/SystemClassicalLight.py.cpp
  ./system/SystemClassicalLight.py.hpp
  ./binding.cpp
)

target_include_directories(bindings PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${Python_INCLUDE_DIRS};${Python_NumPy_INCLUDE_DIRS}"
)

target_link_libraries(bindings
  PUBLIC Eigen3::Eigen
  PRIVATE nanobind-static
)

target_compile_features(bindings PUBLIC cxx_std_20)

set_target_properties(bindings PROPERTIES
  CXX_EXTENSIONS OFF
  POSITION_INDEPENDENT_CODE ON
)

target_compile_options(bindings PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/bigobj /MP>
)

# Assemble objects into libraries
add_library(pairinteraction_backend SHARED
  $<TARGET_OBJECTS:backend>
  $<TARGET_OBJECTS:bindings>
)

target_link_libraries(pairinteraction_backend
  PUBLIC $<$<CXX_COMPILER_ID:MSVC>:${Python_LIBRARIES_GENERALIZED}>
  PRIVATE nanobind-static
)

set_target_properties(pairinteraction_backend PROPERTIES
  PREFIX ""
  INSTALL_RPATH $<$<BOOL:${SKBUILD}>:$<$<PLATFORM_ID:Linux>:$ORIGIN>$<$<PLATFORM_ID:Darwin>:@loader_path>>
)

target_link_options(pairinteraction_backend PRIVATE
  $<$<CXX_COMPILER_ID:AppleClang>:-undefined dynamic_lookup>
  $<$<CXX_COMPILER_ID:Clang,GNU>:-Wl,--unresolved-symbols=ignore-all>
  $<$<CXX_COMPILER_ID:MSVC>:/NODEFAULTLIB:python3${Python_VERSION_MINOR}.lib /DEFAULTLIB:python3.lib /FORCE:UNRESOLVED>
)

# Build an executable for running the tests
if (BUILD_TESTING)
  add_executable(pairinteraction_backend_test
    $<TARGET_OBJECTS:backend>
    main.cpp
  )

  target_link_libraries(pairinteraction_backend_test
    PRIVATE spdlog::spdlog
  )

  # Add tests
  add_test(NAME "pairinteraction_backend_test" COMMAND pairinteraction_backend_test)
endif()
