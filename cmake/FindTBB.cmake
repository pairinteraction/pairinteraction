# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

include(FindPackageHandleStandardArgs)

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
    tbb_config_path = next(p for p in files('tbb-devel') if 'TBBConfig.cmake' in p.name).locate().resolve()
    tbb_library_path = next(p for p in files('tbb') if 'tbb' in p.stem).locate().resolve()
    print(tbb_library_path.parent.parent, tbb_config_path, sep='|')
except PackageNotFoundError:
    sys.exit(1)"
    RESULT_VARIABLE ONEAPI_RESULT
    OUTPUT_VARIABLE ONEAPI_PATHS
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(NOT ONEAPI_RESULT EQUAL 0)
    message(STATUS "Failed to find Intel oneAPI libraries using Python.")
  else()
    string(REPLACE "|" ";" ONEAPI_PATHS_LIST "${ONEAPI_PATHS}")
    list(GET ONEAPI_PATHS_LIST 0 TBB_ROOT)
    list(GET ONEAPI_PATHS_LIST 1 TBB_DIR)
    message(STATUS "TBB root determined to be: ${TBB_ROOT}")
    message(STATUS "TBB package config directory determined to be: ${TBB_DIR}")
    list(APPEND CMAKE_PREFIX_PATH "${TBB_DIR}")
  endif()
else()
  message(STATUS "Python3 interpreter not found; skip discovering Intel oneAPI libraries.")
endif()

find_package(TBB REQUIRED CONFIG)
