# Copyright (c) 2019 Sebastian Weber, Henri Menke. All rights reserved.
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

# Find OpenMP for CXX
#
# OpenMP_CXX_LIBRARY           - the openmp library
# OpenMP_CXX_INCLUDE_DIR       - path including omp.h
# OpenMP::OpenMP_CXX           - target for using OpenMP

include(FindPackageHandleStandardArgs)

find_package(OpenMP QUIET)

if(NOT OpenMP_FOUND AND APPLE)
  find_library(OpenMP_CXX_LIBRARY NAMES omp)
  find_path(OpenMP_CXX_INCLUDE_DIR NAMES omp.h)
  
  if (OpenMP_CXX_LIBRARY AND OpenMP_CXX_INCLUDE_DIR)
    set(OpenMP_CXX_COMPILE_OPTIONS -Xpreprocessor -fopenmp)
    add_library(OpenMP::OpenMP_CXX SHARED IMPORTED)
    set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
      IMPORTED_LOCATION ${OpenMP_CXX_LIBRARY}
      INTERFACE_INCLUDE_DIRECTORIES "${OpenMP_CXX_INCLUDE_DIR}"
      INTERFACE_COMPILE_OPTIONS "${OpenMP_CXX_COMPILE_OPTIONS}"
  )
  endif()
  
  find_package_handle_standard_args(OpenMPCXX DEFAULT_MSG 
    OpenMP_CXX_LIBRARY OpenMP_CXX_INCLUDE_DIR)
    
  mark_as_advanced(OpenMP_CXX_LIBRARY OpenMP_CXX_INCLUDE_DIR)
  
else()
  find_package_handle_standard_args(OpenMPCXX DEFAULT_MSG 
    OpenMP_FOUND)
endif()


