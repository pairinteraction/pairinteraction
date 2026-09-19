# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

# ------------------------------------------------------------
# Generate the backend target
# ------------------------------------------------------------

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
  database/Database.py.hpp
  database/Database.py.cpp
  diagonalize/Diagonalizer.py.cpp
  diagonalize/Diagonalizer.py.hpp
  enums/FloatType.py.hpp
  enums/FloatType.py.cpp
  enums/OperatorType.py.hpp
  enums/OperatorType.py.cpp
  enums/Parity.py.hpp
  enums/Parity.py.cpp
  enums/SorterType.py.hpp
  enums/SorterType.py.cpp
  interfaces/DiagonalizerInterface.py.hpp
  interfaces/DiagonalizerInterface.py.cpp
  interfaces/SorterBuilderInterface.py.hpp
  interfaces/SorterBuilderInterface.py.cpp
  ket/Ket.py.hpp
  ket/Ket.py.cpp
  system/System.py.hpp
  system/System.py.cpp
  tools/run_unit_tests.py.hpp
  tools/run_unit_tests.py.cpp
  bindings.cpp
  LoggerBridge.cpp
  LoggerBridge.hpp
  paths.py.hpp
  paths.py.cpp
  version.py.hpp
  version.py.cpp)

target_link_libraries(backend PRIVATE pairinteraction)

set_target_properties(
  backend
  PROPERTIES OUTPUT_NAME _backend
             INSTALL_RPATH $<$<BOOL:${SKBUILD}>:$<$<PLATFORM_ID:Linux>:$ORIGIN>$<$<PLATFORM_ID:Darwin>:@loader_path>>)

# We mark nanobind's include directory as system include directory to suppress warnings and avoid clang-tidy checks
target_include_directories(backend SYSTEM PRIVATE ${NB_DIR}/include)

# Windows-specific post-build command to copy runtime DLLs
# Windows-specific post-build command to copy runtime DLLs
if(WIN32)
  # OpenSSL's imported targets from vcpkg only define per-config
  # properties for RELEASE and DEBUG (no RELWITHDEBINFO/MINSIZEREL), so
  # the bare IMPORTED_IMPLIB/IMPORTED_LOCATION properties resolve to
  # NOTFOUND for those build types. Look up the config-specific property
  # explicitly, using the Release libs for anything other than Debug.
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(_openssl_config "DEBUG")
  else()
    set(_openssl_config "RELEASE")
  endif()

  get_target_property(OPENSSL_SSL_IMPLIB OpenSSL::SSL IMPORTED_IMPLIB_${_openssl_config})
  get_target_property(OPENSSL_CRYPTO_IMPLIB OpenSSL::Crypto IMPORTED_IMPLIB_${_openssl_config})

  if(NOT OPENSSL_SSL_IMPLIB)
    get_target_property(OPENSSL_SSL_IMPLIB OpenSSL::SSL IMPORTED_LOCATION_${_openssl_config})
  endif()

  if(NOT OPENSSL_SSL_IMPLIB)
    message(FATAL_ERROR "Could not determine OpenSSL::SSL import library for config ${_openssl_config}")
  endif()

  get_filename_component(OPENSSL_LIB_DIR ${OPENSSL_SSL_IMPLIB} DIRECTORY)
  get_filename_component(OPENSSL_BIN_DIR ${OPENSSL_LIB_DIR}/../bin ABSOLUTE)

  file(GLOB OPENSSL_DLLS "${OPENSSL_BIN_DIR}/libssl-*.dll" "${OPENSSL_BIN_DIR}/libcrypto-*.dll")

  if(NOT OPENSSL_DLLS)
    message(FATAL_ERROR "Could not find OpenSSL DLLs in ${OPENSSL_BIN_DIR}")
  endif()

  add_custom_command(
    TARGET backend
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      $<TARGET_RUNTIME_DLLS:backend>
      ${OPENSSL_DLLS}
      $<TARGET_FILE_DIR:backend>
    COMMAND_EXPAND_LISTS)

  install(FILES ${OPENSSL_DLLS} DESTINATION pairinteraction)
  install(FILES $<TARGET_RUNTIME_DLLS:backend> DESTINATION pairinteraction)
endif()

install(IMPORTED_RUNTIME_ARTIFACTS duckdb DESTINATION pairinteraction)

if(WIN32 AND DUCKDB_EXTRA_RUNTIME_DLLS)
  install(FILES ${DUCKDB_EXTRA_RUNTIME_DLLS} DESTINATION pairinteraction)
endif()
install(IMPORTED_RUNTIME_ARTIFACTS TBB::tbb DESTINATION pairinteraction)





# If one cross-compiles, the resulting python package cannot be imported. Thus, a stub file cannot be generated.
if(CMAKE_SYSTEM_PROCESSOR STREQUAL CMAKE_HOST_SYSTEM_PROCESSOR AND NOT CMAKE_CXX_FLAGS MATCHES "sanitize=address")
  # cmake-format: off
  nanobind_add_stub(
    backend_stub
    MODULE _backend
    OUTPUT _backend.pyi
    PYTHON_PATH $<TARGET_FILE_DIR:backend>
    DEPENDS backend
  )
  # cmake-format: on
  # Add mypy ignore comment at the start of the stub file
  install(CODE "file(READ ${CMAKE_CURRENT_BINARY_DIR}/_backend.pyi CONTENT)
          file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/_backend.pyi \"# mypy: ignore-errors\n\n\${CONTENT}\")")
  # Install the .pyi file to the python package directory
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/_backend.pyi DESTINATION pairinteraction)
  # Install the .pyi file also to the python src dir to make static type checking possible with pip editable installs
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/_backend.pyi DESTINATION ${CMAKE_SOURCE_DIR}/src/pairinteraction)
endif()

# ------------------------------------------------------------
# Generate the info target
# ------------------------------------------------------------

# nanobind_add_module info (this internally calls cmake add_library(info ...))
nanobind_add_module(
  info
  STABLE_ABI
  NB_STATIC
  # cmake-format: off
  NB_DOMAIN pairinteraction
  # cmake-format: on
  bindings_info.cpp
  Info.py.hpp
  Info.py.cpp)

target_include_directories(info PUBLIC $<TARGET_PROPERTY:pairinteraction,INTERFACE_INCLUDE_DIRECTORIES>)
target_compile_definitions(info PUBLIC $<TARGET_PROPERTY:pairinteraction,INTERFACE_COMPILE_DEFINITIONS>)

set_target_properties(
  info
  PROPERTIES OUTPUT_NAME _info)

# We mark nanobind's include directory as system include directory to suppress warnings and avoid clang-tidy checks
target_include_directories(info SYSTEM PRIVATE ${NB_DIR}/include)

# If one cross-compiles, the resulting python package cannot be imported. Thus, a stub file cannot be generated.
if(CMAKE_SYSTEM_PROCESSOR STREQUAL CMAKE_HOST_SYSTEM_PROCESSOR AND NOT CMAKE_CXX_FLAGS MATCHES "sanitize=address")
  # cmake-format: off
  nanobind_add_stub(
    info_stub
    MODULE _info
    OUTPUT _info.pyi
    PYTHON_PATH $<TARGET_FILE_DIR:info>
    DEPENDS info
  )
  # cmake-format: on
  # Add mypy ignore comment at the start of the stub file
  install(CODE "file(READ ${CMAKE_CURRENT_BINARY_DIR}/_info.pyi CONTENT)
          file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/_info.pyi \"# mypy: ignore-errors\n\n\${CONTENT}\")")
  # Install the .pyi file to the python package directory
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/_info.pyi DESTINATION pairinteraction)
  # Install the .pyi file also to the python src dir to make static type checking possible with pip editable installs
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/_info.pyi DESTINATION ${CMAKE_SOURCE_DIR}/src/pairinteraction)
endif()

# ------------------------------------------------------------
# Ensure that clang-tidy is really not run on nanobind
# ------------------------------------------------------------

# HACK: Under some circumstances, the marking nanobind's include directory as system include directory does not seem
# to avoid clang-tidy checks. We therefore disable clang-tidy manually for the nanobind-static target that is built
# within nanobind_add_module.
if((Python_VERSION VERSION_LESS 3.12)
  OR (NOT Python_INTERPRETER_ID STREQUAL "Python")
  OR (NOT TARGET Python::SABIModule))
 set_target_properties(nanobind-static PROPERTIES CXX_CLANG_TIDY "")
else()
 set_target_properties(nanobind-static-abi3 PROPERTIES CXX_CLANG_TIDY "")
endif()

# ------------------------------------------------------------
# Install the targets
# ------------------------------------------------------------

install(TARGETS backend info DESTINATION pairinteraction)
