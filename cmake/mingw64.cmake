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

# This defines the relevant dependencies for the MinGW64 toolchain.
#
# Usage: cmake -DCMAKE_TOOLCHAIN_FILE=/path/to/mingw64.cmake

set( CMAKE_SYSTEM_NAME "Windows" )
set( MINGW_PREFIX "x86_64-w64-mingw32" )
set( MINGW_ARCH "x86_64" )

# set compilers
set( CMAKE_C_COMPILER       ${MINGW_PREFIX}-gcc      )
set( CMAKE_CXX_COMPILER     ${MINGW_PREFIX}-g++      )
set( CMAKE_Fortran_COMPILER ${MINGW_PREFIX}-gfortran )

# Set package search directory
set( CMAKE_FIND_ROOT_PATH
      "/usr/${MINGW_PREFIX}"
      "/usr/${MINGW_PREFIX}/sys-root/mingw" )
set( CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY  )
set( CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY  )
set( CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER )
