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

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(pairinteraction CXX)
set(CMAKE_CXX_STANDARD 11)
enable_testing()

# get version info from git

find_package(Git REQUIRED)
execute_process(
  COMMAND "${GIT_EXECUTABLE}" describe --tags --always
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  OUTPUT_VARIABLE VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE )
string(REGEX REPLACE "v(.+)" "\\1" VERSION_WITHOUT_LEADING_V ${VERSION})

# default build type is "Release"

IF(NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Build options

option(WITH_BACKEND  "Build with C++ backend"               ON)
option(WITH_GUI      "Build with Python GUI"                ON)
option(WITH_DATABASE "Generate the quantum defect database" ON)
option(WITH_DOCS     "Generate documentation"               ON)
option(WITH_DMG      "Generate a DMG file (Mac OS X only)"  OFF)

# Mac OS X specific build instructions

set(CMAKE_MACOSX_RPATH 0)
set(CMAKE_MACOSX_GOODIES_PATH "${CMAKE_SOURCE_DIR}/apple")
if( WITH_DMG )
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install")
endif( )

# Build instructions

if( WITH_BACKEND )
  add_subdirectory(calc)
endif( )

if( WITH_DATABASE )
  add_subdirectory(calc/databases)
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
  set(CPACK_PACKAGE_FILE_NAME "pairinteraction-install-linux")
  set(CPACK_PACKAGE_VERSION "${VERSION_WITHOUT_LEADING_V}")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Program for calculation Rydberg interaction potentials")

  # DEB
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "\
    libsqlite3-0, libboost-filesystem-dev, libboost-program-options-dev, \
    libboost-system-dev, libc6, libgcc1, libgsl0ldbl|libgsl2, python3-pint, python3-pyqt5, \
    python3-psutil, python3-pyqt5.qtsvg, python3-scipy, python3-numpy")
  set(CPACK_DEBIAN_PACKAGE_MAINTAINER "PI5")
  set(CPACK_DEBIAN_PACKAGE_SECTION "science")
  set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")

  # RPM
  set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
  set(CPACK_RPM_PACKAGE_REQUIRES "\
    libsqlite3-0, libboost_filesystem1_54_0, libboost_program_options1_54_0, \
    libboost_system1_54_0, glibc, libgcc_s1-gcc49, libgsl0, python3-Pint, python3-qt5, \
    python3-psutil, python3-scipy, python3-numpy")
  set(CPACK_RPM_PACKAGE_VENDOR "PI5")
  set(CPACK_RPM_PACKAGE_GROUP "Productivity/Scientific/Physics")
  set(CPACK_RPM_PACKAGE_PRIORITY "optional")

endif()

include(CPack)

# Print found packages

include(FeatureSummary)
feature_summary(WHAT ALL)
