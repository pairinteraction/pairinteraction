cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

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

else(WIN32)

  set(CPACK_GENERATOR "TGZ;RPM;DEB")

  # General
  set(CPACK_PACKAGE_NAME "pairinteraction")
  set(CPACK_PACKAGE_VERSION "0.0")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Program for calculation Rydberg interaction potentials")

  # DEB
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libsqlite3-0, openmpi-bin, libjsoncpp0, libboost-filesystem-dev, libboost-system-dev, libc6, libgcc1, python3-pint, python3-pyqt5, python3-psutil, python3-pyqt5.qtsvg, python3-scipy, python3-numpy")
  # additional build dependencies: libsqlite3-dev, libeigen3-dev, libjsoncspp-dev, pyqt5-dev-tools, sqlite3, gfortran, libopenmpi-dev
  set(CPACK_DEBIAN_PACKAGE_MAINTAINER "PI5")
  set(CPACK_DEBIAN_PACKAGE_SECTION "science")
  set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
  set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${CMAKE_CURRENT_SOURCE_DIR}/postinst")

  # RPM
  set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
  set(CPACK_RPM_PACKAGE_REQUIRES "libsqlite3-0, openmpi, libjsoncpp1, libboost_filesystem1_54_0, libboost_system1_54_0, glibc, libgcc_s1-gcc49, python3-Pint, python3-qt5, python3-psutil")
  set(CPACK_RPM_PACKAGE_VENDOR "PI5")
  set(CPACK_RPM_PACKAGE_GROUP "Productivity/Scientific/Physics")
  set(CPACK_RPM_PACKAGE_PRIORITY "optional")
  set(CPACK_RPM_POST_INSTALL_SCRIPT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/postinst.rpm")

  include(CPack)

endif(WIN32)
