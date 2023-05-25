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

cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
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

# set version of settings and cache

set(VERSION_SETTINGS "v18")
set(VERSION_CACHE "v18")
set(VERSION_SOFTWARE_DEFAULT "v0.0.0")

# get software version info from git

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

# default build type is "Release with debug info"

IF(NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the build type" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Build options

option(WITH_BACKEND  "Build with C++ backend"               ON)
option(WITH_PYTHON   "Build with SWIG Python interface"     ON)
option(WITH_JULIA    "Build a Julia compatible .so"         OFF)
option(WITH_GUI      "Build with Python GUI"                ON)
option(WITH_DOCS     "Generate documentation"               OFF)
option(WITH_DMG      "Generate a DMG file (Mac OS X only)"  OFF)
option(WITH_COVERAGE "Generate code coverage report"        OFF)
option(WITH_LTO      "Build with link-time optimization"    OFF)
option(WITH_GSL      "Use the GNU scientific library for Whittaker functions" ON)
option(WITH_LAPACK   "Use BLAS and LAPACK to speed up linear algebra" ON)
if(CMAKE_VERSION VERSION_GREATER 3.5.2 AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  option(WITH_CLANG_TIDY "Run Clang-Tidy during compilation" OFF)
endif()

# Enable link-time optimization if available
if(WITH_LTO)
  include(CheckIPOSupported OPTIONAL RESULT_VARIABLE CHECK_IPO_SUPPORTED)
  if(CHECK_IPO_SUPPORTED)
    check_ipo_supported(RESULT IPO_SUPPORTED)
    if(IPO_SUPPORTED)
      set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
  endif()
endif()

# Mac OS X specific build instructions

set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_MACOSX_GOODIES_PATH "${CMAKE_SOURCE_DIR}/apple")
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

# Build instructions

if( WITH_BACKEND )
  add_subdirectory(pairinteraction)
  add_subdirectory(pairinteraction/databases)
  if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    add_subdirectory(testsuite)
  endif()
  if (WITH_JULIA)
    add_subdirectory(pairinteraction.jl/src)
    add_subdirectory(pairinteraction.jl/test)
  endif()
endif( )

if( WITH_GUI )
  add_subdirectory(pairinteraction_gui)
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

configure_file(setup.py.in ${CMAKE_BINARY_DIR}/setup.py @ONLY)
configure_file(setup.cfg.in ${CMAKE_BINARY_DIR}/setup.cfg @ONLY)

if(WIN32)

  # NSIS installer is built externally
  add_subdirectory(win32)

elseif(APPLE AND WITH_DMG)

  find_package(Python3 COMPONENTS Interpreter REQUIRED)

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
    COMMAND ${Python3_EXECUTABLE} "${CMAKE_MACOSX_GOODIES_PATH}/licenseDMG.py" "${PROJECT_BINARY_DIR}/${CPACK_PACKAGE_FILE_NAME}.dmg" "${LICENSES}"
  )

elseif(UNIX)

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
