cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(pairinteraction CXX)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/modules)
set(CMAKE_CXX_COMPILER mpicxx)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

add_definitions(-Wall -Wextra -pedantic -Ofast -march=native -pipe)

# copy databases
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/databases/quantum_defects.db DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

# data type used for scalars
set (SCALAR_DTYPE "complex_t") # "real_t", "complex_t"

# configure a header file to pass some of the CMake settings to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )

# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")

# add the source tree to the search path for include files
include_directories("${PROJECT_SOURCE_DIR}")

# add the executable
file(GLOB pairinteraction_SRCS *.cpp)
add_executable(pairinteraction ${pairinteraction_SRCS})

# add external libraries
find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
else(EIGEN3_FOUND)
  message(FATAL_ERROR "=================\n"
                      "Eigen3 not found!\n"
                      "=================\n")
endif(EIGEN3_FOUND)

find_package(Sqlite3 REQUIRED)
if(SQLITE3_FOUND)
  include_directories(${SQLITE3_INCLUDE_DIR})
  target_link_libraries(pairinteraction ${SQLITE3_LIBRARY})
else(SQLITE3_FOUND)
  message(FATAL_ERROR "==================\n"
                      "Sqlite3 not found!\n"
                      "==================\n")
endif(SQLITE3_FOUND)

find_package(GSL REQUIRED)
if(GSL_FOUND)
  target_link_libraries(pairinteraction ${GSL_LIBRARIES})
else(GSL_FOUND)
  message(FATAL_ERROR "==============\n"
                      "GSL not found!\n"
                      "==============\n")
endif(GSL_FOUND)

find_package(MPI REQUIRED)
if(MPI_FOUND)
  target_link_libraries(pairinteraction ${MPI_LIBRARIES})
else(MPI_FOUND)
  message(FATAL_ERROR "==============\n"
                      "MPI not found!\n"
                      "==============\n")
endif(MPI_FOUND)

find_package(JSONCPP REQUIRED)
if(JSONCPP_FOUND)
  include_directories(${JSONCPP_INCLUDE_DIRS})
  target_link_libraries(pairinteraction ${JSONCPP_LIBRARIES})
else(JSONCPP_FOUND)
  message(FATAL_ERROR "==================\n"
                      "JSONCPP not found!\n"
                      "==================\n")
endif(JSONCPP_FOUND)

# set c++11
set_property(TARGET pairinteraction PROPERTY CXX_STANDARD 11)
set_property(TARGET pairinteraction PROPERTY CXX_STANDARD_REQUIRED ON)
