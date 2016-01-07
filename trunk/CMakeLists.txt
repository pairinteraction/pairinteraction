cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(pairinteraction CXX)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/modules)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_BUILD_TYPE MATCHES Release)
  add_definitions(-Wall -Wextra -pedantic -Ofast -march=native -pipe)
endif(CMAKE_BUILD_TYPE MATCHES Release)

# create output directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)

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
  include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
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
  include_directories(${MPI_CXX_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH})
  target_link_libraries(pairinteraction ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
  set(CMAKE_CXX_FLAGS        ${CMAKE_CXX_FLAGS}        ${MPI_CXX_COMPILE_FLAGS})
  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS})
else(MPI_FOUND)
  message(FATAL_ERROR "==============\n"
                      "MPI not found!\n"
                      "==============\n")
endif(MPI_FOUND)

find_package(Jsoncpp REQUIRED)
if(JSONCPP_FOUND)
  include_directories(${JSONCPP_INCLUDE_DIRS})
  target_link_libraries(pairinteraction ${JSONCPP_LIBRARIES})
else(JSONCPP_FOUND)
  message(FATAL_ERROR "==================\n"
                      "JSONCPP not found!\n"
                      "==================\n")
endif(JSONCPP_FOUND)

find_package(Boost COMPONENTS filesystem system REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(pairinteraction ${Boost_LIBRARIES})
else(Boost_FOUND)
  message(FATAL_ERROR "================\n"
                      "Boost not found!\n"
                      "================\n")
endif(Boost_FOUND)

# set c++11
set_property(TARGET pairinteraction PROPERTY CXX_STANDARD 11)
set_property(TARGET pairinteraction PROPERTY CXX_STANDARD_REQUIRED ON)

# CPack
if(WIN32)

  # TODO: build NSIS installer

else(WIN32)

  install(TARGETS pairinteraction RUNTIME DESTINATION bin)

  set(CPACK_GENERATOR "TGZ;RPM;DEB")

  # General
  set(CPACK_PACKAGE_NAME "pairinteraction")
  set(CPACK_PACKAGE_VERSION "0.0")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Program for calculation Rydberg interaction potentials")

  # DEB
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libsqlite3-0, libgsl0ldbl, openmpi-bin, libjsoncpp0, libboost-filesystem-dev, libboost-system-dev, libc6, libgcc1")
  set(CPACK_DEBIAN_PACKAGE_MAINTAINER "PI5")
  set(CPACK_DEBIAN_PACKAGE_SECTION "science")
  set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")

  # RPM
  set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
  set(CPACK_RPM_PACKAGE_REQUIRES "libsqlite3-0, libgsl0, openmpi, libjsoncpp1, libboost_filesystem1_54_0, libboost_system1_54_0, glibc, libgcc_s1-gcc5")
  set(CPACK_RPM_PACKAGE_VENDOR "PI5")
  set(CPACK_RPM_PACKAGE_GROUP "Productivity/Scientific/Physics")
  set(CPACK_RPM_PACKAGE_PRIORITY "optional")

  include(CPack)

endif(WIN32)
