# SPDX-FileCopyrightText: 2024 Sebastian Weber, Henri Menke, and contributors. All rights reserved.
# SPDX-License-Identifier: GPL-3.0-or-later OR LGPL-3.0-or-later

# Add bindings object
add_library(
  bindings OBJECT
  basis/Basis.py.hpp
  basis/Basis.py.cpp
  database/Database.py.hpp
  database/Database.py.cpp
  diagonalizer/Diagonalizer.py.cpp
  diagonalizer/Diagonalizer.py.hpp
  enums/OperatorType.py.hpp
  enums/TransformationType.py.hpp
  enums/TransformationType.py.cpp
  enums/OperatorType.py.cpp
  interfaces/DiagonalizerInterface.py.hpp
  interfaces/DiagonalizerInterface.py.cpp
  interfaces/TransformationBuilderInterface.py.hpp
  interfaces/TransformationBuilderInterface.py.cpp
  ket/Ket.py.hpp
  ket/Ket.py.cpp
  operator/Operator.py.hpp
  operator/Operator.py.cpp
  system/System.py.hpp
  system/System.py.cpp
  test.py.hpp
  binding.cpp
  test.py.cpp)

target_include_directories(
  bindings
  PRIVATE "$<TARGET_PROPERTY:pairinteraction,INCLUDE_DIRECTORIES>;${Python_INCLUDE_DIRS};${Python_NumPy_INCLUDE_DIRS}")
target_link_libraries(bindings PRIVATE $<$<BOOL:${WITH_COVERAGE}>:coverage::gcov> nanobind-static-pairinteraction)
target_compile_definitions(bindings PRIVATE NB_DOMAIN=pairinteraction)
target_compile_options(bindings PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/bigobj /MP>)

# Assemble objects into backend library
add_library(backend SHARED)

target_link_libraries(
  backend
  PUBLIC $<$<CXX_COMPILER_ID:MSVC>:${Python_LIBRARIES_GENERALIZED}>
  PRIVATE pairinteraction bindings)

set_target_properties(
  backend
  PROPERTIES PREFIX ""
             SUFFIX "${NB_SUFFIX_EXT}"
             INSTALL_RPATH $<$<BOOL:${SKBUILD}>:$<$<PLATFORM_ID:Linux>:$ORIGIN>$<$<PLATFORM_ID:Darwin>:@loader_path>>)

target_link_options(
  backend
  PRIVATE
  $<$<CXX_COMPILER_ID:AppleClang>:-undefined
  dynamic_lookup>
  $<$<CXX_COMPILER_ID:Clang,GNU>:-Wl,--unresolved-symbols=ignore-all>
  $<$<CXX_COMPILER_ID:MSVC>:/NODEFAULTLIB:python3${Python_VERSION_MINOR}.lib
  /DEFAULTLIB:python3.lib
  /FORCE:UNRESOLVED>)

if(WIN32)
  add_custom_command(
    TARGET backend
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:backend> $<TARGET_FILE_DIR:backend>
    COMMAND_EXPAND_LISTS)
endif()

install(IMPORTED_RUNTIME_ARTIFACTS duckdb DESTINATION pairinteraction)
install(TARGETS backend DESTINATION pairinteraction)
if(WIN32)
  install(FILES $<TARGET_RUNTIME_DLLS:backend> DESTINATION pairinteraction)
endif()
