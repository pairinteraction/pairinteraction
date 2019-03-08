# Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
#
# This file is part of the pairinteraction library.
#
# The pairinteraction library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The pairinteraction library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.

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
