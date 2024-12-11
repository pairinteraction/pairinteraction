# SPDX-FileCopyrightText: 2024 Sebastian Weber, Henri Menke, and contributors. All rights reserved.
# SPDX-License-Identifier: GPL-3.0-or-later OR LGPL-3.0-or-later

# nanobind_add_module backend (this internally calls cmake add_library(backend ...))
nanobind_add_module(
  backend
  STABLE_ABI
  NB_STATIC
  # cmake-format: off
  NB_DOMAIN pairinteraction
  # cmake-format: on
  basis/Basis.py.hpp
  basis/Basis.py.cpp
  convenience/matrix_elements.py.hpp
  convenience/matrix_elements.py.cpp
  database/Database.py.hpp
  database/Database.py.cpp
  diagonalizer/Diagonalizer.py.cpp
  diagonalizer/Diagonalizer.py.hpp
  enums/OperatorType.py.hpp
  enums/OperatorType.py.cpp
  enums/Parity.py.hpp
  enums/Parity.py.cpp
  enums/TransformationType.py.hpp
  enums/TransformationType.py.cpp
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
  binding.cpp
  run_unit_tests.py.hpp
  run_unit_tests.py.cpp)

target_link_libraries(backend PRIVATE pairinteraction $<$<BOOL:${WITH_COVERAGE}>:coverage::gcov>)

set_target_properties(
  backend
  PROPERTIES OUTPUT_NAME _backend
             CXX_CLANG_TIDY ""
             INSTALL_RPATH $<$<BOOL:${SKBUILD}>:$<$<PLATFORM_ID:Linux>:$ORIGIN>$<$<PLATFORM_ID:Darwin>:@loader_path>>)

# Windows-specific post-build command to copy runtime DLLs
if(WIN32)
  add_custom_command(
    TARGET backend
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:backend> $<TARGET_FILE_DIR:backend>
    COMMAND_EXPAND_LISTS)
endif()

install(IMPORTED_RUNTIME_ARTIFACTS duckdb DESTINATION pairinteraction/backend)
install(IMPORTED_RUNTIME_ARTIFACTS TBB::tbb DESTINATION pairinteraction/backend)
install(TARGETS backend DESTINATION pairinteraction/backend)

# cmake-format: off
nanobind_add_stub(
  backend_stub
  MODULE _backend
  OUTPUT _backend.pyi
  PYTHON_PATH $<TARGET_FILE_DIR:backend>
  DEPENDS backend
)
# cmake-format: on
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/_backend.pyi DESTINATION pairinteraction/backend)
# install the .pyi file also to the python src dir to make static type checking possible with pip editable installs
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/_backend.pyi DESTINATION ${CMAKE_SOURCE_DIR}/src/pairinteraction/backend)
