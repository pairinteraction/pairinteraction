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

# Find MSMPI
#
# In an effort to have some kind of compatibility with the standard FindMPI we provide
#
# MPI_C_INCLUDE_PATH     - path to header files
# MPI_CXX_INCLUDE_PATH   - ibid.
# MPI_C_LIBRARIES        - path to libraries (architecture dependent)
# MPI_CXX_LIBRARIES      - ibid.

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include (FindPackageHandleStandardArgs)

# Download the MSMPI installer from the Microsoft homepage using curl
# and unpack it using 7-zip.  This is quite ugly and not portable.
set(MS_MPI_SDK_DL "https://download.microsoft.com/download/D/7/B/D7BBA00F-71B7-436B-80BC-4D22F2EE9862/msmpisdk.msi")
set(MSMPI_INCLUDE_PATH ${CMAKE_BINARY_DIR}/MS_MPI)
execute_process( COMMAND curl -C - -L -o ${CMAKE_BINARY_DIR}/msmpisdk.msi ${MS_MPI_SDK_DL} )
execute_process( COMMAND 7z x -aoa -o${MSMPI_INCLUDE_PATH} ${CMAKE_BINARY_DIR}/msmpisdk.msi )

# The library filename is platform dependent. Brrr!
if ("${CMAKE_CXX_COMPILER}" MATCHES "i686" OR "${CMAKE_CXX_COMPILER_ARG1}" MATCHES "i686")
  set( MSMPI_LIBRARIES  "${MSMPI_INCLUDE_PATH}/msmpi.lib" )
elseif ("${CMAKE_CXX_COMPILER}" MATCHES "x86_64" OR "${CMAKE_CXX_COMPILER_ARG1}" MATCHES "x86_64")
  set( MSMPI_LIBRARIES  "${MSMPI_INCLUDE_PATH}/msmpi64.lib" )
endif()

# Make this somewhat compatible to the standard FindMPI
set( MPI_C_INCLUDE_PATH   ${MSMPI_INCLUDE_PATH} )
set( MPI_CXX_INCLUDE_PATH ${MSMPI_INCLUDE_PATH} )

set( MPI_C_LIBRARIES      ${MSMPI_LIBRARIES} )
set( MPI_CXX_LIBRARIES    ${MSMPI_LIBRARIES} )

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args(
  MSMPI DEFAULT_MSG
  MSMPI_INCLUDE_PATH
  MSMPI_LIBRARIES
  )

find_package_handle_standard_args(
  MPI DEFAULT_MSG
  MPI_C_INCLUDE_PATH MPI_CXX_INCLUDE_PATH
  MPI_C_LIBRARIES MPI_CXX_LIBRARIES
  )

mark_as_advanced( 
  MSMPI_INCLUDE_PATH MSMPI_LIBRARIES
  MPI_C_INCLUDE_PATH MPI_CXX_INCLUDE_PATH
  MPI_C_LIBRARIES MPI_CXX_LIBRARIES
  )
