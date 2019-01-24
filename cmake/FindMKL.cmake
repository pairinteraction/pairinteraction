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
        /opt/intel/mkl/lib)

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDE_DIR MKL_LIBRARY)

if (MKL_FOUND)
  if(NOT MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
  endif()
endif()

mark_as_advanced(MKL_LIBRARY MKL_INCLUDE_DIR)
