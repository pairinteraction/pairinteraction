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

if( CMAKE_SYSTEM_NAME MATCHES "Windows" )

  find_program( MAKENSIS NAMES makensis )
  if( MAKENSIS )
    message( STATUS "Found makensis: ${MAKENSIS}" )
  else( )
    message( FATAL_ERROR "Could not find makensis" )
  endif( )

  set( MSMPI_SETUP_DL "https://download.microsoft.com/download/E/8/A/E8A080AF-040D-43FF-97B4-065D4F220301/MSMpiSetup.exe" )
  if( "${CMAKE_CXX_COMPILER}" MATCHES "i686" OR "${CMAKE_CXX_COMPILER_ARG1}" MATCHES "i686" )
    set( MINICONDA_DL "https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86.exe" )
  elseif( "${CMAKE_CXX_COMPILER}" MATCHES "x86_64" OR "${CMAKE_CXX_COMPILER_ARG1}" MATCHES "x86_64" )
    set( MINICONDA_DL "https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe" )
  endif( )

  add_custom_target(win32
    DEPENDS pairinteraction-real pairinteraction-complex database gui
    COMMAND curl -C - -O ${MINICONDA_DL}
    COMMAND curl -C - -O ${MSMPI_SETUP_DL}
    COMMAND ${MAKENSIS} -DARCH=${MINGW_ARCH} ${CMAKE_CURRENT_SOURCE_DIR}/pairinteraction.nsi
    )

endif( )
