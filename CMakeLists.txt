cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(pairinteraction CXX)

set(CMAKE_MACOSX_RPATH 0)

set(CMAKE_COLOR_MAKEFILE ON)

# create output directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)

add_subdirectory(gui2)
add_subdirectory(trunk)

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
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libsqlite3-0, libgsl0ldbl, openmpi-bin, libjsoncpp0, libboost-filesystem-dev, libboost-system-dev, libc6, libgcc1")
  set(CPACK_DEBIAN_PACKAGE_MAINTAINER "PI5")
  set(CPACK_DEBIAN_PACKAGE_SECTION "science")
  set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")

  # RPM
  set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
  set(CPACK_RPM_PACKAGE_REQUIRES "libsqlite3-0, libgsl0, openmpi, libjsoncpp1, libboost_filesystem1_54_0, libboost_system1_54_0, glibc, libgcc_s1-gcc5")
  set(CPACK_RPM_PACKAGE_VENDOR "PI5")
  set(CPACK_RPM_PACKAGE_GROUP "Productivity/Scientific/Physics")
  set(CPACK_RPM_PACKAGE_PRIORITY "optional")

  include(CPack)

endif(WIN32)
