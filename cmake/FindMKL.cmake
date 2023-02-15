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

# Find MKL
#
# MKL_LIBRARIES         - the mkl libraries
# MKL_INCLUDE_DIR       - path including mkl.h

# See https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include(FindPackageHandleStandardArgs)

find_path(MKL_INCLUDE_DIR
  NAMES mkl.h
  PATHS $ENV{MKLROOT}/include
        /opt/intel/mkl/include)

find_library(MKL_LIBRARY
  NAMES mkl_rt
  PATHS $ENV{MKLROOT}/lib
        $ENV{MKLROOT}/lib/intel64
        $ENV{INTEL}/mkl/lib/intel64
        /opt/intel/mkl/lib
        /opt/intel/mkl/lib/intel64)

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDE_DIR MKL_LIBRARY)

if (MKL_FOUND)
  if(NOT MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
  endif()
endif()

mark_as_advanced(MKL_LIBRARY MKL_INCLUDE_DIR)
