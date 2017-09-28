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

# Set compiler options for release builds

if(CMAKE_BUILD_TYPE MATCHES Release)
  if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ox")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -pipe")
  endif()
endif()

# Add current directory to search path

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add the executables

file(GLOB pairinteraction_SRCS *.cpp)
set(wignerSymbols_SRCS wignerSymbols/src/wignerSymbols-cpp.cpp)

add_library(pireal SHARED ${pairinteraction_SRCS} ${wignerSymbols_SRCS})
add_library(picomplex SHARED ${pairinteraction_SRCS} ${wignerSymbols_SRCS})
set_target_properties(picomplex PROPERTIES COMPILE_DEFINITIONS "USE_COMPLEX")

if(MSVC)
  set_source_files_properties(${pairinteraction_SRCS} PROPERTIES COMPILE_FLAGS "/W4 /DNOMINMAX")
  # TODO set_source_files_properties(${wignerSymbols_SRCS} PROPERTIES COMPILE_FLAGS "/wdSOMENUMBER") maybe /wd0169 or /wd0649
else()
  set_source_files_properties(${pairinteraction_SRCS} PROPERTIES COMPILE_FLAGS "-Wall -Wextra -pedantic")
  set_source_files_properties(${wignerSymbols_SRCS} PROPERTIES COMPILE_FLAGS "-Wno-unused-parameter")
endif()

if( WITH_COVERAGE )
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the build type" FORCE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage -fprofile-arcs -ftest-coverage")
  link_libraries(gcov)
endif( )

# Use OpenMP

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Add external libraries

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

find_package(Sqlite3 REQUIRED)
include_directories(${SQLITE3_INCLUDE_DIR})
list(APPEND LIBRARIES ${SQLITE3_LIBRARY})

find_package(Boost COMPONENTS filesystem system program_options serialization REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
list(APPEND LIBRARIES ${Boost_LIBRARIES})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
list(APPEND LIBRARIES ${GSL_LIBRARIES})

find_package(ZeroMQ REQUIRED)
include_directories(${ZEROMQ_INCLUDE_DIR})
list(APPEND LIBRARIES ${ZEROMQ_LIBRARY})

target_link_libraries(pireal ${LIBRARIES})
target_link_libraries(picomplex ${LIBRARIES})

# Generate interface with SWIG

set( DATATYPE "real" )
configure_file( Interface.i.cmakein ${CMAKE_CURRENT_BINARY_DIR}/interface_real.i )
set( DATATYPE "complex" )
configure_file( Interface.i.cmakein ${CMAKE_CURRENT_BINARY_DIR}/interface_complex.i )

find_package(SWIG REQUIRED)
find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)
find_package(Numpy REQUIRED)

include(${SWIG_USE_FILE})
include_directories(${PYTHON_INCLUDE_PATH})
include_directories(SYSTEM ${NUMPY_INCLUDE_DIR})

set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/interface_real.i PROPERTIES CPLUSPLUS ON SWIG_FLAGS "-py3")
swig_add_module(pireal python ${CMAKE_CURRENT_BINARY_DIR}/interface_real.i)
swig_link_libraries(pireal pireal ${PYTHON_LIBRARIES})

set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/interface_complex.i PROPERTIES CPLUSPLUS ON SWIG_FLAGS "-py3;-DUSE_COMPLEX")
swig_add_module(picomplex python ${CMAKE_CURRENT_BINARY_DIR}/interface_complex.i)
swig_link_libraries(picomplex picomplex ${PYTHON_LIBRARIES})
set_source_files_properties( ${swig_generated_file_fullname} PROPERTIES COMPILE_DEFINITIONS "USE_COMPLEX")

# Include test cases

add_subdirectory( unit_test )

# CPack

if( WITH_DMG )
  find_package(PythonInterp 3 REQUIRED)

  install(TARGETS _pireal LIBRARY DESTINATION pairinteraction)
  install(TARGETS _picomplex LIBRARY DESTINATION pairinteraction)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pireal.py DESTINATION pairinteraction)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/picomplex.py DESTINATION pairinteraction)
  install(TARGETS pireal LIBRARY DESTINATION pairinteraction)
  install(TARGETS picomplex LIBRARY DESTINATION pairinteraction)
  
  set(bin1 \${CMAKE_INSTALL_PREFIX}/pairinteraction/_pireal.so)
  set(bin2 \${CMAKE_INSTALL_PREFIX}/pairinteraction/_picomplex.so)
  set(bin3 \${CMAKE_INSTALL_PREFIX}/pairinteraction/libpireal.dylib)
  set(bin4 \${CMAKE_INSTALL_PREFIX}/pairinteraction/libpicomplex.dylib)
  
  install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_MACOSX_GOODIES_PATH}/standalone.py ${bin1} ${bin2} ${bin3} ${bin4})")

elseif ( NOT WIN32 )

  install(TARGETS _pireal LIBRARY DESTINATION share/pairinteraction/pairinteraction)
  install(TARGETS _picomplex LIBRARY DESTINATION share/pairinteraction/pairinteraction)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pireal.py DESTINATION share/pairinteraction/pairinteraction)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/picomplex.py DESTINATION share/pairinteraction/pairinteraction)
  install(TARGETS pireal LIBRARY DESTINATION share/pairinteraction/pairinteraction)
  install(TARGETS picomplex LIBRARY DESTINATION share/pairinteraction/pairinteraction)

endif( )
