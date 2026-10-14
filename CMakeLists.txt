# SPDX-FileCopyrightText: 2024 PairInteraction Developers
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
  $ pip install . --group tests --group docs
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

include(FetchContent)

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

# Our vcpkg triplet (cmake/triplets/x64-windows.cmake) builds only the release versions of the dependencies. Because
# objects that are compiled against the debug runtime are binary incompatible with them, we use the release runtime for
# all build types, including Debug. For the same reason, imported targets that provide debug versions (e.g., TBB::tbb from
# tbb-devel) must use their release versions in Debug builds, too. We detect a release-only triplet by the missing debug
# directory instead of by the triplet name, because vcpkg's built-in triplet of the same name builds debug versions.
# Other setups keep CMake's defaults.
if(MSVC
   AND DEFINED VCPKG_TARGET_TRIPLET
   AND IS_DIRECTORY "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin"
   AND NOT EXISTS "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/debug"
   AND NOT DEFINED CMAKE_MSVC_RUNTIME_LIBRARY)
  set(CMAKE_MSVC_RUNTIME_LIBRARY MultiThreadedDLL)
  if(NOT DEFINED CMAKE_MAP_IMPORTED_CONFIG_DEBUG)
    # If none of the listed configurations is provided, CMake treats the imported target as not found. Thus, we list
    # all other configurations, too: NoConfig (packages exported without a build type), the empty entry (targets without
    # configuration-specific properties), and Debug last (for packages that only provide debug versions).
    set(CMAKE_MAP_IMPORTED_CONFIG_DEBUG "Release;RelWithDebInfo;MinSizeRel;NoConfig;;Debug")
  endif()
  # In Debug builds, vcpkg's toolchain copies the DLLs of the dependencies from the non-existent debug/bin directory, i.e.,
  # not at all. $<TARGET_RUNTIME_DLLS> covers only DLLs of imported shared targets but not, e.g., the ICU DLLs needed by
  # duckdb. Thus, the post-build steps that copy the runtime DLLs additionally copy the release DLLs in Debug builds.
  file(GLOB VCPKG_RELEASE_DLLS_FOR_DEBUG "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin/*.dll")
endif()

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

# Find duckdb
set(DUCKDB_EXPECTED_VERSION 1.4.4) # keep in sync with pyproject.toml
find_package(DuckDB QUIET)
# DuckDB reports its version with a "v" prefix (e.g., "v1.4.4"), which CMake's version comparison treats as older than
# any "1.x.y". Thus, we cannot pass the version to find_package and compare the versions ourselves.
if(DuckDB_FOUND)
  string(REGEX REPLACE "^v" "" DUCKDB_FOUND_VERSION "${DuckDB_VERSION}")
  if(NOT DUCKDB_FOUND_VERSION VERSION_EQUAL DUCKDB_EXPECTED_VERSION)
    message(WARNING "Found DuckDB ${DuckDB_VERSION} in ${DuckDB_DIR}, but version ${DUCKDB_EXPECTED_VERSION} is "
                    "expected to match the duckdb Python package.")
  endif()
endif()
if(NOT DuckDB_FOUND AND (APPLE OR UNIX))
  if(APPLE)
    set(DUCKDB_URL https://install.duckdb.org/v${DUCKDB_EXPECTED_VERSION}/libduckdb-osx-universal.zip)
    set(DUCKDB_LIBNAME libduckdb.dylib)
    set(DUCKDB_SHA256 5e6d79f5fb86bbd4f24d59cee3bc38f113d1433761df35337e3c01c62eeafe26)
  elseif(UNIX)
    # You can manually compile DuckDB and use it by calling CMake with
    # -DFETCHCONTENT_SOURCE_DIR_DUCKDB=/path/to/folder/containing/libduckdb.so/and/duckdb.hpp.
    set(DUCKDB_URL https://install.duckdb.org/v${DUCKDB_EXPECTED_VERSION}/libduckdb-linux-amd64.zip)
    set(DUCKDB_LIBNAME libduckdb.so)
    set(DUCKDB_SHA256 09cc288295964d897b47665d1898e16e8ef176cae9ea615797fc136eae15bd5d)
  endif()
  FetchContent_Declare(
    duckdb
    URL ${DUCKDB_URL} CONFIGURE_COMMAND "" BUILD_COMMAND ""
    URL_HASH SHA256=${DUCKDB_SHA256})
  FetchContent_GetProperties(duckdb)
  if(NOT duckdb_POPULATED)
    FetchContent_Populate(duckdb)
  endif()
  add_library(duckdb SHARED IMPORTED)
  target_include_directories(duckdb INTERFACE ${duckdb_SOURCE_DIR})
  set_target_properties(
    duckdb PROPERTIES IMPORTED_LOCATION ${duckdb_SOURCE_DIR}/${DUCKDB_LIBNAME}
                      INTERFACE_SYSTEM_INCLUDE_DIRECTORIES $<TARGET_PROPERTY:duckdb,INTERFACE_INCLUDE_DIRECTORIES>)
elseif(NOT DuckDB_FOUND)
    message(FATAL_ERROR "DuckDB not found and automatic fallback is only supported on macOS/Linux.")
endif()

# Doctest
# We use the dev version because it works with recent cmake versions that forbid a version < 3.5 for cmake_minimum_required
FetchContent_Declare(
  doctest
  GIT_REPOSITORY https://github.com/doctest/doctest
  GIT_TAG 3a01ec37828affe4c9650004edb5b304fb9d5b75 # v0.20.0
)
# Add the EXCLUDE_FROM_ALL option to avoid adding files of doctest to the wheel. In the future, with cmake 3.28, we
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
set(HTTPLIB_USE_ZSTD_IF_AVAILABLE OFF)
set(HTTPLIB_USE_NON_BLOCKING_GETADDRINFO OFF)
FetchContent_Declare(
  httplib
  GIT_REPOSITORY https://github.com/yhirose/cpp-httplib
  GIT_TAG bd95e67c234930cd6d6bb11309588c5462c63cec # v0.30.1
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

list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")

include(CTest)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  add_subdirectory(src/cpp/tests)
endif()

# Print found packages
include(FeatureSummary)
feature_summary(WHAT ALL)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
