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

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(DMG "FALSE" CACHE BOOL "indicates whether a .dmg file should be created")

if(DMG)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install")
endif(DMG)

project(pairinteraction CXX)

set(CMAKE_MACOSX_RPATH 0)

set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_VERBOSE_MAKEFILE OFF)

IF(NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Release ... FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_subdirectory(gui)
add_subdirectory(calc)

# CPack
if(WIN32)

  # TODO: build NSIS installer

elseif(APPLE AND DMG)

  # http://stackoverflow.com/questions/13420117/how-to-get-a-path-of-the-app-bundle-in-a-cmake-script-to-supply-it-to-the-fixup
  # https://github.com/mgalloy/cpack-cheatsheets/blob/master/cpack-variables.rst
  # http://stackoverflow.com/questions/12302836/renaming-cpack-automatic-target
  # http://stackoverflow.com/questions/24163778/howto-add-custom-target-that-depends-on-make-install
  # https://cmake.org/pipermail/cmake/2015-September/061533.html
  
  set(CPACK_GENERATOR "Bundle")
  
  set(CPACK_PACKAGE_FILE_NAME "pairinteraction-install-osx")
  set(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/apple/pairinteraction.icns")
  set(CPACK_BUNDLE_NAME "pairinteraction")
  set(CPACK_BUNDLE_ICON "${CMAKE_SOURCE_DIR}/apple/pairinteraction.icns")
  set(CPACK_BUNDLE_PLIST "${CMAKE_SOURCE_DIR}/apple/Info.plist")
  set(CPACK_BUNDLE_STARTUP_COMMAND "${CMAKE_SOURCE_DIR}/apple/launcher.sh")
  
  include(CPack)

elseif(UNIX)

  set(CPACK_GENERATOR "TGZ;RPM;DEB")

  # General
  set(CPACK_PACKAGE_NAME "pairinteraction")
  set(CPACK_PACKAGE_VERSION "0.0")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Program for calculation Rydberg interaction potentials")

  # DEB
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libsqlite3-0, openmpi-bin, libboost-filesystem-dev, libboost-system-dev, libc6, libgcc1, python3-pint, python3-pyqt5, python3-psutil, python3-pyqt5.qtsvg, python3-scipy, python3-numpy")
  # additional build dependencies: libsqlite3-dev, libeigen3-dev, pyqt5-dev-tools, sqlite3, gfortran, libopenmpi-dev
  set(CPACK_DEBIAN_PACKAGE_MAINTAINER "PI5")
  set(CPACK_DEBIAN_PACKAGE_SECTION "science")
  set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
  set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${CMAKE_CURRENT_SOURCE_DIR}/postinst")

  # RPM
  set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
  set(CPACK_RPM_PACKAGE_REQUIRES "libsqlite3-0, openmpi, libboost_filesystem1_54_0, libboost_system1_54_0, glibc, libgcc_s1-gcc49, python3-Pint, python3-qt5, python3-psutil")
  set(CPACK_RPM_PACKAGE_VENDOR "PI5")
  set(CPACK_RPM_PACKAGE_GROUP "Productivity/Scientific/Physics")
  set(CPACK_RPM_PACKAGE_PRIORITY "optional")
  set(CPACK_RPM_POST_INSTALL_SCRIPT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/postinst.rpm")

  include(CPack)

endif()
