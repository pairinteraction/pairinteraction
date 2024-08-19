# Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
#
# This file is part of the pairinteraction library.
#
# The pairinteraction library is free software: you can redistribute it and/or modify it under the terms of the GNU
# Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# The pairinteraction library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along with the pairinteraction library. If
# not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
message(STATUS "CMake version: ${CMAKE_VERSION}")

project(pairinteraction CXX)
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
  $ pip install .[gui,doc]
=====================================================================
Where the arguments 'gui' (installs the graphical user interface) \
and 'doc' (installs dependencies for building the documentation) \
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
option(WITH_DOC "Generate documentation" OFF)
option(WITH_DMG "Generate a DMG file (Mac OS X only)" OFF) # TODO treat building a windows installer and dmg the same
option(WITH_COVERAGE "Generate code coverage report" OFF)

# Windows specific build instructions
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# Mac OS X specific build instructions
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_MACOSX_GOODIES_PATH "${CMAKE_SOURCE_DIR}/extra/apple")
if(WITH_DMG)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install")
endif()

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
find_package(doctest REQUIRED)
find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(TBB REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(
  Python REQUIRED
  COMPONENTS Interpreter Development.Module NumPy
  OPTIONAL_COMPONENTS Development.SABIModule)

find_package(httplib REQUIRED)
if(NOT httplib_OpenSSL_FOUND)
  find_package(OpenSSL REQUIRED COMPONENTS SSL Crypto)
  add_compile_definitions(CPPHTTPLIB_OPENSSL_SUPPORT)
endif()

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
if(MSVC)
  string(REPLACE "python3${Python_VERSION_MINOR}.lib" "python3.lib" Python_LIBRARIES_GENERALIZED ${Python_LIBRARIES})
  message(STATUS "Unpatched Python libraries: ${Python_LIBRARIES}")
  message(STATUS "  Patched Python libraries: ${Python_LIBRARIES_GENERALIZED}")
endif()
include(FetchContent)
FetchContent_Declare(
  nanobind
  GIT_REPOSITORY https://github.com/wjakob/nanobind.git
  GIT_TAG 8d7f1ee0621c17fa370b704b2100ffa0243d5bfb # v2.0.0
)
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
  set(DUCKDB_URL https://github.com/duckdb/duckdb/releases/download/v1.0.0/libduckdb-windows-amd64.zip)
  set(DUCKDB_LIBNAME duckdb.dll)
  set(DUCKDB_MD5 86585ecd2429a9c47cb7fcb9d4b3a168)
elseif(APPLE)
  set(DUCKDB_URL https://github.com/duckdb/duckdb/releases/download/v1.0.0/libduckdb-osx-universal.zip)
  set(DUCKDB_LIBNAME libduckdb.dylib)
  set(DUCKDB_MD5 70cd2d330ef721c5854cdad51fceb8bc)
elseif(UNIX)
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
add_library(duckdb::duckdb ALIAS duckdb)

# Add subdirectories for building the project
add_subdirectory(src/cpp)

add_subdirectory(src/cpp/bindings)

include(CTest)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  add_subdirectory(src/cpp/tests)
endif()

if(WITH_DOC)
  add_subdirectory(doc)
endif()

# Concatenate licenses
file(READ ${CMAKE_SOURCE_DIR}/LICENSE-LGPL.txt LICENSE_LGPL)
file(READ ${CMAKE_SOURCE_DIR}/LICENSE-GPL.txt LICENSE_GPL)
file(READ ${CMAKE_SOURCE_DIR}/licenses-3rd-party/LICENSE-MKL.txt LICENSE_MKL)

string(
  CONCAT LICENSE_ALL ${LICENSE_LGPL}
         "\n\n\n=============================================================================\n\n\n" ${LICENSE_GPL}
         "\n\n\n=============================================================================\n\n\n" ${LICENSE_MKL})

file(WRITE ${PROJECT_BINARY_DIR}/LICENSE-ALL.txt ${LICENSE_ALL})

# Make packages for different platforms

if(APPLE AND WITH_DMG)

  find_package(
    Python
    COMPONENTS Interpreter
    REQUIRED)

  install(FILES ${CMAKE_SOURCE_DIR}/LICENSE-LGPL.txt ${CMAKE_SOURCE_DIR}/LICENSE-GPL.txt
                ${CMAKE_SOURCE_DIR}/licenses-3rd-party/LICENSE-MKL.txt DESTINATION licenses)

  set(CPACK_GENERATOR "Bundle")

  set(CPACK_PACKAGE_FILE_NAME "pairinteraction-install-osx")
  set(CPACK_PACKAGE_ICON "${CMAKE_MACOSX_GOODIES_PATH}/pairinteraction.icns")
  set(CPACK_BUNDLE_NAME "pairinteraction")
  set(CPACK_BUNDLE_ICON "${CMAKE_MACOSX_GOODIES_PATH}/pairinteraction.icns")
  set(CPACK_BUNDLE_PLIST "${CMAKE_MACOSX_GOODIES_PATH}/Info.plist")
  set(CPACK_BUNDLE_STARTUP_COMMAND "${CMAKE_MACOSX_GOODIES_PATH}/launcher.sh")

  add_custom_target(
    license COMMAND ${Python_EXECUTABLE} ${CMAKE_MACOSX_GOODIES_PATH}/licenseDMG.py
                    ${PROJECT_BINARY_DIR}/${CPACK_PACKAGE_FILE_NAME}.dmg ${PROJECT_BINARY_DIR}/LICENSE-ALL.txt)

elseif(UNIX AND NOT SKBUILD)

  add_subdirectory(extra/flatpak)

  install(FILES ${CMAKE_SOURCE_DIR}/README.md DESTINATION share/doc/${CMAKE_PROJECT_NAME}/)
  install(
    FILES ${PROJECT_BINARY_DIR}/LICENSE-ALL.txt
    DESTINATION share/doc/${CMAKE_PROJECT_NAME}/
    RENAME copyright)

  set(CPACK_GENERATOR "TGZ;RPM;DEB")

  # General
  if(NOT CPACK_PACKAGE_FILE_NAME)
    set(CPACK_PACKAGE_FILE_NAME "pairinteraction-install-linux")
  endif()
  string(REGEX REPLACE "\\.(deb|rpm|tar\\.gz)$" "" CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_FILE_NAME})
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
  set(CPACK_RPM_PACKAGE_REQUIRES "python3-Pint, python3-qt5, python3-scipy, python3-numpy, python3-pyqtgraph")
  set(CPACK_RPM_PACKAGE_VENDOR "PI5")
  set(CPACK_RPM_PACKAGE_GROUP "Productivity/Scientific/Physics")
  set(CPACK_RPM_PACKAGE_PRIORITY "optional")

endif()

include(CPack)

# Print found packages
include(FeatureSummary)
feature_summary(WHAT ALL)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# HACK: Copy the python directories to the build directory # TODO re
if(NOT SKBUILD)
  file(COPY ${CMAKE_SOURCE_DIR}/src/pairinteraction DESTINATION ${CMAKE_BINARY_DIR})
endif()
