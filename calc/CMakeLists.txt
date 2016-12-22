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

if(CMAKE_BUILD_TYPE MATCHES Release)
  add_definitions(-Wall -Wextra -pedantic -Ofast -pipe)
endif(CMAKE_BUILD_TYPE MATCHES Release)

# Add the executables

file(GLOB pairinteraction_SRCS *.cpp)
set(pairinteraction_matrixelement_SRCS ${pairinteraction_SRCS})
set(wignerSymbols_SRCS wignerSymbols/src/wignerSymbols-cpp.cpp)
list(REMOVE_ITEM pairinteraction_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/mainMatrixelement.cpp)
list(REMOVE_ITEM pairinteraction_matrixelement_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

add_executable(pairinteraction-real ${pairinteraction_SRCS} ${wignerSymbols_SRCS})
add_executable(pairinteraction-complex ${pairinteraction_SRCS} ${wignerSymbols_SRCS})
add_executable(matrixelements ${pairinteraction_matrixelement_SRCS} ${wignerSymbols_SRCS})
set_target_properties(pairinteraction-complex PROPERTIES COMPILE_DEFINITIONS "USE_COMPLEX")

set_source_files_properties(${wignerSymbols_SRCS} PROPERTIES COMPILE_FLAGS -Wno-unused-parameter)

# Add external libraries

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})


find_package(Sqlite3 REQUIRED)
include_directories(${SQLITE3_INCLUDE_DIR})
target_link_libraries(pairinteraction-real ${SQLITE3_LIBRARY})
target_link_libraries(pairinteraction-complex ${SQLITE3_LIBRARY})
target_link_libraries(matrixelements ${SQLITE3_LIBRARY})


if( WIN32 )
  find_package(MSMPI REQUIRED)
  set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc ${CMAKE_EXE_LINKER_FLAGS}")
else( )
  find_package(MPI REQUIRED)
endif( )
include_directories(${MPI_CXX_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH})
target_link_libraries(pairinteraction-real ${MPI_CXX_LIBRARIES})
target_link_libraries(pairinteraction-complex ${MPI_CXX_LIBRARIES})
target_link_libraries(matrixelements ${MPI_CXX_LIBRARIES})
set(CMAKE_CXX_FLAGS        ${CMAKE_CXX_FLAGS}        ${MPI_CXX_COMPILE_FLAGS})
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS})


find_package(Boost COMPONENTS filesystem system program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(pairinteraction-real ${Boost_LIBRARIES})
target_link_libraries(pairinteraction-complex ${Boost_LIBRARIES})
target_link_libraries(matrixelements ${Boost_LIBRARIES})


find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
target_link_libraries(pairinteraction-real ${GSL_LIBRARIES})
target_link_libraries(pairinteraction-complex ${GSL_LIBRARIES})
target_link_libraries(matrixelements ${GSL_LIBRARIES})

# Set c++11

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(FATAL "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# CPack

if( WITH_DMG )

  install(TARGETS pairinteraction-real RUNTIME DESTINATION .)
  install(TARGETS pairinteraction-complex RUNTIME DESTINATION .)
  install(TARGETS matrixelements RUNTIME DESTINATION .)
  
  set(bin1 \${CMAKE_INSTALL_PREFIX}/pairinteraction-real)
  set(bin2 \${CMAKE_INSTALL_PREFIX}/pairinteraction-complex)
  set(bin3 \${CMAKE_INSTALL_PREFIX}/matrixelements)
  
  install(CODE "execute_process(COMMAND ${PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/apple/standalone.py ${bin1} ${bin2} ${bin3})")

else( )

  install(TARGETS pairinteraction-real RUNTIME DESTINATION share/pairinteraction/calc)
  install(TARGETS pairinteraction-complex RUNTIME DESTINATION share/pairinteraction/calc)
  install(TARGETS matrixelements RUNTIME DESTINATION share/pairinteraction/calc)

endif( )
