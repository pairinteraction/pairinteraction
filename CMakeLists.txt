# Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

project(pairinteraction CXX)
set(CMAKE_CXX_STANDARD 11)
enable_testing()

# Enable link-time optimization if available
if(POLICY CMP0069)
  cmake_policy(SET CMP0069 NEW)
  include(CheckIPOSupported OPTIONAL RESULT_VARIABLE CHECK_IPO_SUPPORTED)
  if(CHECK_IPO_SUPPORTED)
    check_ipo_supported(RESULT IPO_SUPPORTED)
    if(IPO_SUPPORTED)
      set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
  endif()
endif()

set(LIBNAME "libpairinteraction")

# set version of settings and cache

set(VERSION_SETTINGS "v17")
set(VERSION_CACHE "v16")

# get software version info from git

find_package(Git)

if(GIT_FOUND)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" describe --tags --always
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE VERSION_SOFTWARE
    OUTPUT_STRIP_TRAILING_WHITESPACE )
else()
  message(WARNING "Because Git was not found, the version of the pairinteraction software cannot be determined.")
  set(VERSION_SOFTWARE "v0.0.0")
endif()

string(REGEX REPLACE "v(.+)" "\\1" VERSION_WITHOUT_LEADING_V ${VERSION_SOFTWARE})

# default build type is "Release"

IF(NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Build options

option(WITH_BACKEND  "Build with C++ backend"               ON)
option(WITH_PYTHON   "Build with SWIG Python interface"     ON)
option(WITH_GUI      "Build with Python GUI"                ON)
option(WITH_DOCS     "Generate documentation"               ON)
option(WITH_DMG      "Generate a DMG file (Mac OS X only)"  OFF)
option(WITH_COVERAGE "Generate code coverage report"        OFF)
if(CMAKE_VERSION VERSION_GREATER 3.5.2 AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  option(WITH_CLANG_TIDY "Run Clang-Tidy during compilation" OFF)
endif()

# Mac OS X specific build instructions

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_MACOSX_GOODIES_PATH "${CMAKE_SOURCE_DIR}/apple")
if( WITH_DMG )
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
    set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}"
      "-checks=-*,bugprone-*,misc-*,modernize-*,-modernize-use-using,performance-*,readability-*,cert-*,-cert-flp30-c,clang-analyzer-*,-clang-analyzer-security.FloatLoopCounter"
      "-header-filter=.*"
      "-warnings-as-errors=*")
  endif()
endif()

# Build instructions

if( WITH_BACKEND )
  add_subdirectory(${LIBNAME})
  add_subdirectory(${LIBNAME}/databases)
  if(WITH_PYTHON)
    add_subdirectory(testsuite)
  endif()
endif( )

if( WITH_GUI )
  add_subdirectory(gui)
endif( )

if ( WITH_DOCS )
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

if(WIN32)

  # NSIS installer is built externally
  add_subdirectory(win32)

elseif(APPLE AND WITH_DMG)

  find_package(PythonInterp 3 REQUIRED)

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
    COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_MACOSX_GOODIES_PATH}/licenseDMG.py" "${PROJECT_BINARY_DIR}/${CPACK_PACKAGE_FILE_NAME}.dmg" "${LICENSES}"
  )
  
elseif(UNIX)

  install(
    FILES ${CMAKE_SOURCE_DIR}/LICENSE.txt
    DESTINATION share/doc/${CMAKE_PROJECT_NAME}/
    COMPONENT applications
    RENAME copyright
  )

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
    "python3-pint, python3-pyqt5, python3-psutil, python3-pyqt5.qtsvg, python3-scipy, python3-numpy")
  set(CPACK_DEBIAN_PACKAGE_MAINTAINER "PI5")
  set(CPACK_DEBIAN_PACKAGE_SECTION "science")
  set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")

  # RPM
  set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
  set(CPACK_RPM_PACKAGE_AUTOREQPROV ON)
  set(CPACK_RPM_PACKAGE_REQUIRES
    "python3-Pint, python3-qt5, python3-psutil, python3-scipy, python3-numpy")
  set(CPACK_RPM_PACKAGE_VENDOR "PI5")
  set(CPACK_RPM_PACKAGE_GROUP "Productivity/Scientific/Physics")
  set(CPACK_RPM_PACKAGE_PRIORITY "optional")

endif()

include(CPack)

# Print found packages

include(FeatureSummary)
feature_summary(WHAT ALL)
