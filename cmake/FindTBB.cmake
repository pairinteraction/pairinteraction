# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

include(FindPackageHandleStandardArgs)

# Remember how this module was called, the nested find_package calls below overwrite these variables
set(TBB_IS_REQUIRED "${TBB_FIND_REQUIRED}")
set(TBB_IS_QUIET "${TBB_FIND_QUIETLY}")

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
    message(STATUS "Failed to find the 'tbb-devel' Python package using ${Python3_EXECUTABLE}.")
  else()
    string(REPLACE "|" ";" ONEAPI_PATHS_LIST "${ONEAPI_PATHS}")
    list(GET ONEAPI_PATHS_LIST 0 TBB_ROOT)
    list(GET ONEAPI_PATHS_LIST 1 TBB_CONFIG_FILE)
    # TBB_DIR must be the directory that contains the config file, not the config file itself
    get_filename_component(TBB_DIR "${TBB_CONFIG_FILE}" DIRECTORY)
    message(STATUS "TBB root determined to be: ${TBB_ROOT}")
    message(STATUS "TBB package config directory determined to be: ${TBB_DIR}")
    list(APPEND CMAKE_PREFIX_PATH "${TBB_DIR}")
  endif()
else()
  message(STATUS "Python3 interpreter not found; skip discovering Intel oneAPI libraries.")
endif()

find_package(TBB QUIET CONFIG)

set(TBB_FIND_REQUIRED "${TBB_IS_REQUIRED}")
set(TBB_FIND_QUIETLY "${TBB_IS_QUIET}")

find_package_handle_standard_args(
  TBB
  REQUIRED_VARS TBB_DIR
  VERSION_VAR TBB_VERSION
  REASON_FAILURE_MESSAGE
    "TBB is obtained from the 'tbb-devel' Python package. Install the build requirements into the Python environment \
that CMake uses (${Python3_EXECUTABLE}) by running 'pip install -r .build_requirements.txt'.")
