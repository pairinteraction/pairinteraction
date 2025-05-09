# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
message(STATUS "CMake version: ${CMAKE_VERSION}")

project(pairinteraction CXX C)
if(NOT SKBUILD)
  message(
    WARNING
      "\
This CMake file is mainly meant to be executed using 'scikit-build' \
or by persons developing the C++ backend. Running it directly might \
not produce the desired result. If you are a user trying to install \
this package, consider using the following install command inside \
the root directory of the repository:
=====================================================================
  $ pip install .[tests,docs]
=====================================================================
Where the arguments 'tests' (installs dependencies for running tests) \
and 'docs' (installs dependencies for building the documentation) \
are optional and can be omitted. \

Optionally: Set up a virtual environment and activate it to install \
the package in an isolated environment. We recommend using 'uv' \
(https://pypi.org/project/uv) for this like so:
=====================================================================
  $ uv venv .venv
  $ . .venv/bin/activate
=====================================================================
Now you can use `uv pip install .`. \

For more details and advanced options, please refer to installation \
instructions in the documentation.")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Default build type is "Release with debug info"
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      RelWithDebInfo
      CACHE STRING "Choose the build type" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Build options
option(WITH_COVERAGE "Generate code coverage report" OFF)

# Windows specific build instructions
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# Mac OS X specific build instructions
set(CMAKE_MACOSX_RPATH TRUE)

# General build instructions
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Target for code coverage
if(WITH_COVERAGE)
  if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING "-DWITH_COVERAGE=ON implies -DCMAKE_BUILD_TYPE=Debug")
  endif()
  set(CMAKE_BUILD_TYPE
      "Debug"
      CACHE STRING "Choose the build type" FORCE)
  add_library(coverage::gcov INTERFACE IMPORTED)
  target_compile_options(coverage::gcov INTERFACE $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:--coverage -fprofile-arcs
                                                  -ftest-coverage>)
  target_link_libraries(coverage::gcov INTERFACE $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:gcov>)
endif()

# Find dependencies
find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(fmt REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(
  Python REQUIRED
  COMPONENTS Interpreter Development.Module NumPy
  OPTIONAL_COMPONENTS Development.SABIModule)
find_package(TBB REQUIRED)

# Try to find MKL
set(CMAKE_C_COMPILER_ID "${CMAKE_CXX_COMPILER_ID}") # Needed by MKLConfig to find gnu_thread
set(MKL_ARCH
    "intel64"
    CACHE STRING "MKL architecture")
set(MKL_LINK
    "dynamic"
    CACHE STRING "MKL linkage")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(MKL_THREADING
      "sequential"
      CACHE STRING "MKL threading backend")
else()
  set(MKL_THREADING
      "tbb_thread"
      CACHE STRING "MKL threading backend")
endif()
set(MKL_INTERFACE
    "lp64"
    CACHE STRING "MKL interface")
find_package(MKL)

# Try to find lapacke if MKL is not available
if(NOT MKL_FOUND)
  find_package(LAPACKE)
endif()

# Find nanobind
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE NB_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

# Find duckdb  # keep the version used here in sync with pyproject.toml
include(FetchContent)
if(WIN32)
  set(DUCKDB_URL https://github.com/duckdb/duckdb/releases/download/v1.1.3/libduckdb-windows-amd64.zip)
  set(DUCKDB_LIBNAME duckdb.dll)
  set(DUCKDB_MD5 422a6317ad652f995930828d457a715a)
elseif(APPLE)
  set(DUCKDB_URL https://github.com/duckdb/duckdb/releases/download/v1.1.3/libduckdb-osx-universal.zip)
  set(DUCKDB_LIBNAME libduckdb.dylib)
  set(DUCKDB_MD5 1430f9f3427a245f2d6d89ccd1a5196c)
elseif(UNIX)
  # Note that the latest binary release of DuckDB does not use the CXX11 ABI because it still uses gcc4 (see also
  # https://stackoverflow.com/questions/77938452/duckdb-c-linking-error-undefined-reference and
  # https://github.com/duckdb/duckdb/blob/51529635f58039f0caadce9bd3304ed53c6b0b5d/.github/workflows/LinuxRelease.yml).
  # This is why we use v0.8.1 here. We recommend to manually compile a more recent version of DuckDB and use it by
  # calling CMake with -DFETCHCONTENT_SOURCE_DIR_DUCKDB=/path/to/folder/containing/libduckdb.so/and/duckdb.hpp.
  set(DUCKDB_URL https://github.com/duckdb/duckdb/releases/download/v0.8.1/libduckdb-linux-amd64.zip)
  set(DUCKDB_LIBNAME libduckdb.so)
  set(DUCKDB_MD5 249fad0221153fb99af58d4a150d5fb0)
endif()
FetchContent_Declare(
  duckdb
  URL ${DUCKDB_URL} CONFIGURE_COMMAND "" BUILD_COMMAND ""
  URL_HASH MD5=${DUCKDB_MD5})
FetchContent_GetProperties(duckdb)
if(NOT duckdb_POPULATED)
  FetchContent_Populate(duckdb)
endif()
add_library(duckdb SHARED IMPORTED)
target_include_directories(duckdb INTERFACE ${duckdb_SOURCE_DIR})
set_target_properties(
  duckdb PROPERTIES IMPORTED_LOCATION ${duckdb_SOURCE_DIR}/${DUCKDB_LIBNAME}
                    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES $<TARGET_PROPERTY:duckdb,INTERFACE_INCLUDE_DIRECTORIES>)
if(WIN32)
  set_target_properties(duckdb PROPERTIES IMPORTED_IMPLIB ${duckdb_SOURCE_DIR}/duckdb.lib)
endif()

# Doctest
# We use the dev version because it works with recent cmake versions that forbid a version < 3.5 for cmake_minimum_required
FetchContent_Declare(
  doctest
  GIT_REPOSITORY https://github.com/doctest/doctest
  GIT_TAG 3a01ec37828affe4c9650004edb5b304fb9d5b75 # v0.20.0
)
# Add the EXCLUDE_FROM_ALL option to avoid adding files of httplib to the wheel. In the future, with cmake 3.28, we
# can add this argument directly to FetchContent_Declare.
if(NOT doctest_POPULATED)
  FetchContent_Populate(doctest)
  add_subdirectory(${doctest_SOURCE_DIR} ${doctest_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
FetchContent_MakeAvailable(doctest)

# Httplib
set(HTTPLIB_REQUIRE_OPENSSL ON)
set(HTTPLIB_USE_BROTLI_IF_AVAILABLE OFF)
set(HTTPLIB_USE_ZLIB_IF_AVAILABLE OFF)
FetchContent_Declare(
  httplib
  GIT_REPOSITORY https://github.com/yhirose/cpp-httplib
  GIT_TAG 787a34ad7f01f20922a237d5142aae469828be72 # v0.20.0
)
# Add the EXCLUDE_FROM_ALL option to avoid adding files of httplib to the wheel. In the future, with cmake 3.28, we
# can add this argument directly to FetchContent_Declare.
FetchContent_GetProperties(httplib)
if(NOT httplib_POPULATED)
  FetchContent_Populate(httplib)
  add_subdirectory(${httplib_SOURCE_DIR} ${httplib_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
FetchContent_MakeAvailable(httplib)

# Find miniz
FetchContent_Declare(
  miniz
  URL https://github.com/richgel999/miniz/releases/download/3.0.2/miniz-3.0.2.zip
  URL_HASH MD5=0604f14151944ff984444b04c5c760e5
)
FetchContent_GetProperties(miniz)
if(NOT miniz_POPULATED)
  FetchContent_Populate(miniz)
endif()
add_library(miniz STATIC ${miniz_SOURCE_DIR}/miniz.h  ${miniz_SOURCE_DIR}/miniz.c)
if(UNIX)
  target_compile_definitions(miniz PUBLIC _LARGEFILE64_SOURCE)
endif()
target_include_directories(miniz SYSTEM PUBLIC ${miniz_SOURCE_DIR})

# Find Cpptrace
find_package(cpptrace QUIET)
if(NOT cpptrace_FOUND)
  FetchContent_Declare(
    cpptrace
    SYSTEM
    GIT_REPOSITORY https://github.com/jeremy-rifkin/cpptrace.git
    GIT_TAG 6689d14c203eed390ae7bb64f56a983cfd7dff9c # v0.7.5
  )
  # Add the EXCLUDE_FROM_ALL option to avoid adding files of cpptrace to the wheel. In the future, with cmake 3.28, we
  # can add this argument directly to FetchContent_Declare.
  FetchContent_GetProperties(cpptrace)
  if(NOT cpptrace_POPULATED)
    FetchContent_Populate(cpptrace)
    add_subdirectory(${cpptrace_SOURCE_DIR} ${cpptrace_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
  FetchContent_MakeAvailable(cpptrace)
  set_target_properties(cpptrace-lib PROPERTIES CXX_CLANG_TIDY "")
endif()

# Add subdirectories for building the project
add_subdirectory(src/cpp)
add_subdirectory(src/cpp/bindings)
add_subdirectory(src/cpp/docs)

list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")

include(CTest)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  add_subdirectory(src/cpp/tests)
endif()

# Print found packages
include(FeatureSummary)
feature_summary(WHAT ALL)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
