# Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
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

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(basics CXX)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Ofast -pipe")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add executable

add_executable(basics basics.cpp)

# Enable OpenMP if supported

find_package(OpenMP)
if (${OPENMP_FOUND})
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

################################################################################
### TODO create a FindPairinteraction.cmake file ############################### 
################################################################################

# Add external libraries needed by pairinteraction

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../../../cmake) # TODO

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

# Add the pairinteraction library

#target_compile_definitions(basics PRIVATE USE_COMPLEX) # if used, libpicomplex.dylib has to be linked
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../)
list(APPEND LIBRARIES "/Applications/pairinteraction.app/Contents/Resources/libpairinteraction/libpireal.dylib")

################################################################################

# Link boost filesystem

find_package(Boost COMPONENTS filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
list(APPEND LIBRARIES ${Boost_LIBRARIES})

# Link libraries

target_link_libraries(basics ${LIBRARIES})
