# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

include(FindPackageHandleStandardArgs)

if(MKL_THREADING STREQUAL "tbb_thread")
  find_package(TBB REQUIRED)
endif()

find_package(
  Python3
  COMPONENTS Interpreter
  QUIET)
if(Python3_FOUND)
  execute_process(
    COMMAND
      ${Python3_EXECUTABLE} -c "import sys
from importlib.metadata import files, PackageNotFoundError
try:
    mkl_config_path = next(p for p in files('mkl-devel') if 'MKLConfig.cmake' in p.name).locate().resolve()
    mkl_library_path = next(p for p in files('mkl') if 'mkl_core' in p.stem).locate().resolve()
    print(mkl_library_path.parent.parent, mkl_config_path, sep='|')
except PackageNotFoundError:
    sys.exit(1)"
    RESULT_VARIABLE ONEAPI_RESULT
    OUTPUT_VARIABLE ONEAPI_PATHS
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(NOT ONEAPI_RESULT EQUAL 0)
    message(STATUS "Failed to find Intel oneAPI libraries using Python.")
  else()
    string(REPLACE "|" ";" ONEAPI_PATHS_LIST "${ONEAPI_PATHS}")
    list(GET ONEAPI_PATHS_LIST 0 MKL_ROOT)
    list(GET ONEAPI_PATHS_LIST 1 MKL_DIR)
    message(STATUS "MKL root determined to be: ${MKL_ROOT}")
    message(STATUS "MKL package config directory determined to be: ${MKL_DIR}")
    list(APPEND CMAKE_PREFIX_PATH "${MKL_DIR}")
  endif()
else()
  message(STATUS "Python3 interpreter not found; skip discovering Intel oneAPI libraries.")
endif()

find_package(MKL QUIET CONFIG)
if(MKL_CONSIDERED_CONFIGS)
  set(target MKL::MKL)

  find_package_handle_standard_args(MKL CONFIG_MODE)

  if(OMP_LIBRARY)
    get_filename_component(OMP_LIBRARY_PATH "${OMP_LIBRARY}" DIRECTORY)
    get_filename_component(OMP_LIBRARY_PATH "${OMP_LIBRARY_PATH}" REALPATH)
    set_property(
      TARGET MKL::MKL
      APPEND
      PROPERTY LINK_OPTIONS "$<$<CXX_COMPILER_ID:Clang,GNU>:LINKER:-rpath=${OMP_LIBRARY_PATH}>")
  endif()
else()
  find_package(PkgConfig QUIET)
  if(PKG_CONFIG_FOUND)
    if(MKL_LINK STREQUAL "sdl")
      pkg_search_module(MKL QUIET IMPORTED_TARGET mkl-sdl)
    elseif(MKL_THREADING STREQUAL "sequential")
      pkg_search_module(MKL QUIET IMPORTED_TARGET mkl-${MKL_LINK}-${MKL_INTERFACE}-seq)
    elseif(MKL_THREADING STREQUAL "intel_thread")
      pkg_search_module(MKL QUIET IMPORTED_TARGET mkl-${MKL_LINK}-${MKL_INTERFACE}-iomp)
    elseif(MKL_THREADING STREQUAL "gnu_thread")
      pkg_search_module(MKL QUIET IMPORTED_TARGET mkl-${MKL_LINK}-${MKL_INTERFACE}-gomp)
    elseif(MKL_THREADING STREQUAL "tbb_thread")
      pkg_search_module(MKL QUIET IMPORTED_TARGET mkl-${MKL_LINK}-${MKL_INTERFACE}-tbb)
    endif()
    set(target PkgConfig::MKL)
  endif()

  find_package_handle_standard_args(
    MKL
    REQUIRED_VARS MKL_FOUND
    VERSION_VAR MKL_VERSION)

  if(MKL_FOUND)
    add_library(MKL::MKL ALIAS PkgConfig::MKL)
  endif()
endif()

if(MKL_FOUND AND WIN32)
  # https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-windows/2023-1/contents-of-the-redist-intel64-directory.html
  file(
    GLOB
    MKL_DLLS
    # Threading layer
    "${MKL_ROOT}/bin/mkl_${MKL_THREADING}*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_${MKL_THREADING}*.dll"
    # Computational layer
    "${MKL_ROOT}/bin/mkl_core*.dll"
    "${MKL_ROOT}/bin/mkl_def*.dll"
    "${MKL_ROOT}/bin/mkl_mc*.dll"
    "${MKL_ROOT}/bin/mkl_avx*.dll"
    "${MKL_ROOT}/bin/mkl_vml_def*.dll"
    "${MKL_ROOT}/bin/mkl_vml_mc*.dll"
    "${MKL_ROOT}/bin/mkl_vml_avx*.dll"
    "${MKL_ROOT}/bin/mkl_vml_cmpt*.dll"
    "${MKL_ROOT}/bin/libimalloc.dll"
    "${MKL_ROOT}/redist/intel64/mkl_core*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_def*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_mc*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_avx*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_vml_def*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_vml_mc*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_vml_avx*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_vml_cmpt*.dll"
    "${MKL_ROOT}/redist/intel64/libimalloc.dll")

  if(NOT MKL_DLLS)
    message(WARNING "MKL DLLs could not be located! Dynamic lookup will probably fail.")
  endif()

  if(MKL_THREADING STREQUAL "intel_thread")
    # Implementation detail of MKLConfig.cmake, the pkg-config files do not expose these variables
    if(EXISTS "${OMP_DLL_DIR}/${OMP_DLLNAME}")
      list(APPEND MKL_DLLS ${OMP_DLL_DIR}/${OMP_DLLNAME})
    else()
      message(WARNING "Intel OpenMP DLL cloud not be located! Dynamic lookup will probably fail.")
    endif()
  endif()

  # Because there are no suitable implibs for the MKL DLLs, we use the one from mkl_core
  get_target_property(MKL_CORE_IMPORTED_IMPLIB MKL::mkl_core IMPORTED_IMPLIB)

  foreach(FILE ${MKL_DLLS})
    message(STATUS "Including MKL DLL: ${FILE}")
    string(MD5 FILE_HASH "${FILE}")
    add_library("mkl_${FILE_HASH}" SHARED IMPORTED)
    set_target_properties("mkl_${FILE_HASH}" PROPERTIES IMPORTED_LOCATION "${FILE}" IMPORTED_IMPLIB
                                                                                    "${MKL_CORE_IMPORTED_IMPLIB}")
    target_link_libraries(${target} INTERFACE "mkl_${FILE_HASH}")
  endforeach()
endif()
