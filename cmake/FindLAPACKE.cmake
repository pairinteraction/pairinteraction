# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

include(FindPackageHandleStandardArgs)
include(CheckCXXSymbolExists)

find_package(LAPACK QUIET)

if(LAPACK_FOUND)
  # Some systems have LAPACKE in a separate library for some reason. The following attempts to detect that. Static
  # linking is preferred.
  find_library(
    LAPACKE_LIBRARY
    NAMES liblapacke.a liblapacke.dylib lapacke
    PATH_SUFFIXES lib
    PATHS /usr/local/opt/lapack /opt/homebrew/opt/lapack)
  find_path(
    LAPACKE_INCLUDE_DIR
    NAMES lapacke.h
    PATH_SUFFIXES include
    PATHS /usr/include/openblas /usr/local/opt/lapack /opt/homebrew/opt/lapack)

  set(CMAKE_REQUIRED_INCLUDES_SAVED ${CMAKE_REQUIRED_INCLUDES})
  set(CMAKE_REQUIRED_LIBRARIES_SAVED ${CMAKE_REQUIRED_LIBRARIES})
  if(LAPACKE_INCLUDE_DIR)
    set(CMAKE_REQUIRED_INCLUDES "${LAPACKE_INCLUDE_DIR}")
  endif()
  if(LAPACKE_LIBRARY)
    set(CMAKE_REQUIRED_LIBRARIES "${LAPACKE_LIBRARY};${LAPACK_LIBRARIES}")
  else()
    set(CMAKE_REQUIRED_LIBRARIES "${LAPACK_LIBRARIES}")
  endif()
  check_cxx_symbol_exists(LAPACKE_dsyevd "lapacke.h" LAPACKE_WORKS)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES_SAVED})
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_SAVED})

  if(LAPACKE_WORKS)
    if(LAPACKE_LIBRARY)
      set_property(TARGET LAPACK::LAPACK PROPERTY INTERFACE_LINK_LIBRARIES "${LAPACKE_LIBRARY};${LAPACK_LIBRARIES}")
    endif()
    if(LAPACKE_INCLUDE_DIR)
      set_property(
        TARGET LAPACK::LAPACK
        APPEND
        PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${LAPACKE_INCLUDE_DIR}")
    endif()
    add_library(LAPACKE::LAPACKE ALIAS LAPACK::LAPACK)
  endif()
endif()

find_package_handle_standard_args(
  LAPACKE
  VERSION_VAR LAPACK_VERSION
  REQUIRED_VARS LAPACK_FOUND LAPACKE_WORKS
  HANDLE_COMPONENTS)
