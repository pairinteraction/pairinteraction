# Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
#
# This file is part of the pairinteraction library.
#
# The pairinteraction library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The pairinteraction library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
message(STATUS "CMake version: ${CMAKE_VERSION}")

project(pairinteraction CXX)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()

# Work around buggy CMake ignoring CXX_EXTENSIONS OFF
# https://gitlab.kitware.com/cmake/cmake/-/issues/22224
if(POLICY CMP0128)
  cmake_policy(SET CMP0128 NEW)
endif()

# Set version of settings and cache
set(VERSION_SETTINGS "v18")
set(VERSION_CACHE "v18")
set(VERSION_SOFTWARE_DEFAULT "v0.0.0")

# Get software version info from git
if(NOT VERSION_SOFTWARE)
  find_package(Git)
  if(GIT_FOUND)
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" describe --abbrev=0 --tags
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      RESULT_VARIABLE GIT_DESCRIBE_RESULT
      OUTPUT_VARIABLE GIT_DESCRIBE_OUTPUT
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(GIT_DESCRIBE_RESULT EQUAL 0)
      set(VERSION_SOFTWARE "${GIT_DESCRIBE_OUTPUT}" CACHE STRING "Set software version" FORCE)
    else()
      message(WARNING "Git failed to determine the version of the pairinteraction software.")
      set(VERSION_SOFTWARE "${VERSION_SOFTWARE_DEFAULT}" CACHE STRING "Set software version" FORCE)
    endif()
  else()
    message(WARNING "Because Git was not found, the version of the pairinteraction software cannot be determined.")
    set(VERSION_SOFTWARE "${VERSION_SOFTWARE_DEFAULT}" CACHE STRING "Set software version" FORCE)
  endif()
endif()

string(REGEX REPLACE "v(.+)" "\\1" VERSION_WITHOUT_LEADING_V ${VERSION_SOFTWARE})
string(REGEX REPLACE "-(pre)?alpha" "a" VERSION_WITHOUT_LEADING_V_COMPACT_ALPHA ${VERSION_WITHOUT_LEADING_V})

# Default build type is "Release with debug info"
IF(NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the build type" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Build options
option(WITH_BACKEND  "Build with C++ backend"               ON)
option(WITH_PYTHON   "Build with SWIG Python interface"     ON)
option(WITH_GUI      "Build with Python GUI"                ON)
option(WITH_DOC      "Generate documentation"               OFF)
option(WITH_BENCH    "Compile the benchmarks"               OFF)
option(WITH_DMG      "Generate a DMG file (Mac OS X only)"  OFF)
option(WITH_COVERAGE "Generate code coverage report"        OFF)
option(WITH_GSL      "Use the GNU scientific library for Whittaker functions" ON)
option(WITH_LAPACK   "Use BLAS and LAPACK to speed up linear algebra" ON)

# Windows specific build instructions
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# Mac OS X specific build instructions
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_MACOSX_GOODIES_PATH "${CMAKE_SOURCE_DIR}/extra/apple")
if( WITH_DMG )
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install")
endif( )

# Target for the unit tests
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C ${CMAKE_BUILD_TYPE})

# Perform static analysis using clang-tidy
if(WITH_CLANG_TIDY)
  string(REGEX REPLACE "^([1-9]+\\.[0-9]+).*$" "\\1" CLANG_MINOR_VERSION "${CMAKE_CXX_COMPILER_VERSION}")
  find_program(
    CLANG_TIDY_EXE
    NAMES "clang-tidy-${CLANG_MINOR_VERSION}" "clang-tidy"
    DOC "Path to clang-tidy executable"
    )
  if(NOT CLANG_TIDY_EXE)
    message(STATUS "clang-tidy not found.")
  else()
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
  endif()
endif()

# Find dependencies # TODO add all dependencies
find_package(doctest REQUIRED)
find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Python REQUIRED
  COMPONENTS Interpreter Development.Module NumPy
  OPTIONAL_COMPONENTS Development.SABIModule
)

find_package(httplib COMPONENTS OpenSSL)
if (NOT httplib_FOUND)
  find_package(httplib REQUIRED)
  find_package(OpenSSL REQUIRED COMPONENTS SSL Crypto)
  add_compile_definitions(CPPHTTPLIB_OPENSSL_SUPPORT)
endif()

# Find MKL
set(CMAKE_C_COMPILER_ID "${CMAKE_CXX_COMPILER_ID}") # Needed by MKLConfig to find gnu_thread
set(MKL_ARCH "intel64" CACHE STRING "MKL architecture")
set(MKL_LINK "dynamic" CACHE STRING "MKL linkage")
set(MKL_THREADING "sequential" CACHE STRING "MKL threading backend")
set(MKL_INTERFACE "lp64" CACHE STRING "MKL interface")
find_package(MKL)

# Find nanobind
if(MSVC)
  string(REPLACE "python3${Python_VERSION_MINOR}.lib" "python3.lib" Python_LIBRARIES_GENERALIZED ${Python_LIBRARIES})
  message(STATUS "Unpatched Python libraries: ${Python_LIBRARIES}")
  message(STATUS "  Patched Python libraries: ${Python_LIBRARIES_GENERALIZED}")
endif()
include(FetchContent)
FetchContent_Declare(nanobind
  GIT_REPOSITORY https://github.com/wjakob/nanobind.git
  GIT_TAG        72831bb)
FetchContent_GetProperties(nanobind)
if(NOT nanobind_POPULATED)
  FetchContent_Populate(nanobind)
endif()
list(APPEND CMAKE_PREFIX_PATH ${nanobind_SOURCE_DIR})
find_package(nanobind CONFIG REQUIRED)
nanobind_build_library(nanobind-static)
set_target_properties(nanobind-static PROPERTIES CXX_CLANG_TIDY "")

# Find duckdb
if(WIN32)
  set(DUCKDB_URL https://github.com/duckdb/duckdb/releases/download/v0.10.2/libduckdb-windows-amd64.zip)
  set(DUCKDB_LIBNAME duckdb.dll)
elseif(APPLE)
  set(DUCKDB_URL https://github.com/duckdb/duckdb/releases/download/v0.10.2/libduckdb-osx-universal.zip)
  set(DUCKDB_LIBNAME libduckdb.dylib)
elseif(UNIX)
  set(DUCKDB_URL https://github.com/duckdb/duckdb/releases/download/v0.8.1/libduckdb-linux-amd64.zip)
  set(DUCKDB_LIBNAME libduckdb.so)
endif()
FetchContent_Declare(duckdb
    URL ${DUCKDB_URL}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
)
FetchContent_GetProperties(duckdb)
if(NOT duckdb_POPULATED)
  FetchContent_Populate(duckdb)
endif()
add_library(duckdb SHARED IMPORTED)
target_include_directories(duckdb INTERFACE ${duckdb_SOURCE_DIR})
set_target_properties(duckdb PROPERTIES
  IMPORTED_LOCATION ${duckdb_SOURCE_DIR}/${DUCKDB_LIBNAME}
  INTERFACE_SYSTEM_INCLUDE_DIRECTORIES $<TARGET_PROPERTY:duckdb,INTERFACE_INCLUDE_DIRECTORIES>)
if(WIN32)
  set_target_properties(duckdb PROPERTIES IMPORTED_IMPLIB ${duckdb_SOURCE_DIR}/duckdb.lib)
endif()
add_library(duckdb::duckdb ALIAS duckdb)

# Build instructions
if( WITH_BACKEND )
  add_subdirectory(pairinteraction_backend)
  if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    add_subdirectory(test/pairinteraction_backend/cpp)
    add_subdirectory(test/pairinteraction_gui) # TODO use pytest
    add_subdirectory(test/pairinteraction) # TODO use pytest
  endif()
endif( )

if(WITH_BENCH)
endif()

if( WITH_GUI )
  add_subdirectory(pairinteraction_gui)
endif( )

if ( WITH_DOC )
  add_subdirectory(doc)
endif( )

# Concatenate licenses
set(LICENSES "${PROJECT_BINARY_DIR}/LICENSES.txt")

file(READ "${CMAKE_SOURCE_DIR}/LICENSE.txt" LICENSESTXT)
file(READ "${CMAKE_SOURCE_DIR}/LICENSE-3RD-PARTY.txt" 3RDPARTYTXT)

if("${3RDPARTYTXT}" MATCHES "(^.*)\n\n#*\nLICENSE TEXTS")
  string(CONCAT LICENSESTXT "${LICENSESTXT}" "-----------------------------------------------------------------------------\n\n" "${CMAKE_MATCH_1}")
endif()

file(WRITE "${LICENSES}" "${LICENSESTXT}")

# Make packages for different platforms
configure_file(setup.py.in ${CMAKE_BINARY_DIR}/setup.py @ONLY)
configure_file(setup.cfg.in ${CMAKE_BINARY_DIR}/setup.cfg @ONLY)

if(WIN32)

  # NSIS installer is built externally
  add_subdirectory(extra/win32)

elseif(APPLE AND WITH_DMG)

  find_package(Python COMPONENTS Interpreter REQUIRED)

  install(
    FILES ${CMAKE_SOURCE_DIR}/LICENSE.txt ${CMAKE_SOURCE_DIR}/LICENSE-3RD-PARTY.txt
    DESTINATION licenses
  )

  set(CPACK_GENERATOR "Bundle")

  set(CPACK_PACKAGE_FILE_NAME "pairinteraction-install-osx")
  set(CPACK_PACKAGE_ICON "${CMAKE_MACOSX_GOODIES_PATH}/pairinteraction.icns")
  set(CPACK_BUNDLE_NAME "pairinteraction")
  set(CPACK_BUNDLE_ICON "${CMAKE_MACOSX_GOODIES_PATH}/pairinteraction.icns")
  set(CPACK_BUNDLE_PLIST "${CMAKE_MACOSX_GOODIES_PATH}/Info.plist")
  set(CPACK_BUNDLE_STARTUP_COMMAND "${CMAKE_MACOSX_GOODIES_PATH}/launcher.sh")

  add_custom_target(
    license
    COMMAND ${Python_EXECUTABLE} "${CMAKE_MACOSX_GOODIES_PATH}/licenseDMG.py" "${PROJECT_BINARY_DIR}/${CPACK_PACKAGE_FILE_NAME}.dmg" "${LICENSES}"
  )

elseif(UNIX AND NOT SKBUILD)

  add_subdirectory(extra/flatpak)

  install(FILES ${CMAKE_SOURCE_DIR}/README.md
    DESTINATION share/doc/${CMAKE_PROJECT_NAME}/)
  install(FILES ${CMAKE_SOURCE_DIR}/LICENSE.txt
    DESTINATION share/doc/${CMAKE_PROJECT_NAME}/
    RENAME copyright)

  set(CPACK_GENERATOR "TGZ;RPM;DEB")

  # General
  if(NOT CPACK_PACKAGE_FILE_NAME)
    set(CPACK_PACKAGE_FILE_NAME "pairinteraction-install-linux")
  endif()
  string(REGEX REPLACE "\\.(deb|rpm|tar\\.gz)$" "" CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_FILE_NAME})
  set(CPACK_PACKAGE_VERSION "${VERSION_WITHOUT_LEADING_V}")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Program for calculation Rydberg interaction potentials")

  # DEB
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
  set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS
    "python3-pint, python3-pyqt5, python3-pyqt5.qtsvg, python3-scipy, python3-numpy, python3-pyqtgraph")
  set(CPACK_DEBIAN_PACKAGE_MAINTAINER "PI5")
  set(CPACK_DEBIAN_PACKAGE_SECTION "science")
  set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")

  # RPM
  set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
  set(CPACK_RPM_PACKAGE_AUTOREQPROV ON)
  set(CPACK_RPM_PACKAGE_REQUIRES
    "python3-Pint, python3-qt5, python3-scipy, python3-numpy, python3-pyqtgraph")
  set(CPACK_RPM_PACKAGE_VENDOR "PI5")
  set(CPACK_RPM_PACKAGE_GROUP "Productivity/Scientific/Physics")
  set(CPACK_RPM_PACKAGE_PRIORITY "optional")

endif()

include(CPack)

# Print found packages
include(FeatureSummary)
feature_summary(WHAT ALL)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Software version: ${VERSION_SOFTWARE}")
if(WITH_GSL)
  message(STATUS "License compatibility: GPL v3")
else()
  message(STATUS "License compatibility: LGPL v3")
endif()

# HACK: Copy the "pairinteraction" directory to the build directory # TODO re
if(NOT SKBUILD)
  file(COPY ${CMAKE_SOURCE_DIR}/pairinteraction DESTINATION ${CMAKE_BINARY_DIR})
endif()
