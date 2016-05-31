cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(CMAKE_MACOSX_RPATH 0)

add_subdirectory(databases)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/modules)

if(CMAKE_BUILD_TYPE MATCHES Release)
  add_definitions(-Wall -Wextra -pedantic -Ofast -march=native -pipe)
endif(CMAKE_BUILD_TYPE MATCHES Release)

# copy databases
# file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/databases/quantum_defects.db DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

# add the executables
file(GLOB pairinteraction_SRCS *.cpp)
set(pairinteraction_matrixelement_SRCS ${pairinteraction_SRCS})
list(REMOVE_ITEM pairinteraction_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/mainMatrixelement.cpp)
list(REMOVE_ITEM pairinteraction_matrixelement_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

add_executable(pairinteraction-real ${pairinteraction_SRCS})
add_executable(pairinteraction-complex ${pairinteraction_SRCS})
add_executable(matrixelements ${pairinteraction_matrixelement_SRCS})
set_target_properties(pairinteraction-complex PROPERTIES COMPILE_DEFINITIONS "USE_COMPLEX")

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
  target_link_libraries(pairinteraction-real ${SQLITE3_LIBRARY})
  target_link_libraries(pairinteraction-complex ${SQLITE3_LIBRARY})
  target_link_libraries(matrixelements ${SQLITE3_LIBRARY})
else(SQLITE3_FOUND)
  message(FATAL_ERROR "==================\n"
                      "Sqlite3 not found!\n"
                      "==================\n")
endif(SQLITE3_FOUND)

if("${CMAKE_CXX_COMPILER}" MATCHES "mingw32" OR "${CMAKE_CXX_COMPILER_ARG1}" MATCHES "mingw32")

  set(MS_MPI_SDK_DL "https://download.microsoft.com/download/D/7/B/D7BBA00F-71B7-436B-80BC-4D22F2EE9862/msmpisdk.msi")
  set(MS_MPI_INCLUDE_PATH ${CMAKE_BINARY_DIR}/MS_MPI)
  execute_process(COMMAND curl -C - -L -o ${CMAKE_BINARY_DIR}/msmpisdk.msi ${MS_MPI_SDK_DL})
  execute_process(COMMAND 7z x -aoa -o${MS_MPI_INCLUDE_PATH} ${CMAKE_BINARY_DIR}/msmpisdk.msi)
  include_directories(${MS_MPI_INCLUDE_PATH})
  target_link_libraries(pairinteraction-real "${MS_MPI_INCLUDE_PATH}/msmpi.lib")
  target_link_libraries(pairinteraction-complex "${MS_MPI_INCLUDE_PATH}/msmpi.lib")
  target_link_libraries(matrixelements "${MS_MPI_INCLUDE_PATH}/msmpi.lib")

else()

  find_package(MPI REQUIRED)
  if(MPI_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH})
    target_link_libraries(pairinteraction-real ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
    target_link_libraries(pairinteraction-complex ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
    target_link_libraries(matrixelements ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
    set(CMAKE_CXX_FLAGS        ${CMAKE_CXX_FLAGS}        ${MPI_CXX_COMPILE_FLAGS})
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS})
  else(MPI_FOUND)
    message(FATAL_ERROR "==============\n"
                        "MPI not found!\n"
                        "==============\n")
  endif(MPI_FOUND)

endif()

find_package(Jsoncpp REQUIRED)
if(JSONCPP_FOUND)
  include_directories(${JSONCPP_INCLUDE_DIRS})
  target_link_libraries(pairinteraction-real ${JSONCPP_LIBRARIES})
  target_link_libraries(pairinteraction-complex ${JSONCPP_LIBRARIES})
  target_link_libraries(matrixelements ${JSONCPP_LIBRARIES})
else(JSONCPP_FOUND)
  message(FATAL_ERROR "==================\n"
                      "JSONCPP not found!\n"
                      "==================\n")
endif(JSONCPP_FOUND)

find_package(Boost COMPONENTS filesystem system REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(pairinteraction-real ${Boost_LIBRARIES})
  target_link_libraries(pairinteraction-complex ${Boost_LIBRARIES})
  target_link_libraries(matrixelements ${Boost_LIBRARIES})
else(Boost_FOUND)
  message(FATAL_ERROR "================\n"
                      "Boost not found!\n"
                      "================\n")
endif(Boost_FOUND)

# set c++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# add library wignerSymbols
add_subdirectory(wignerSymbols)
include_directories(wignerSymbols/include)
target_link_libraries(pairinteraction-real wignerSymbols)
target_link_libraries(pairinteraction-complex wignerSymbols)
target_link_libraries(matrixelements wignerSymbols)


# CPack

install(TARGETS pairinteraction-real RUNTIME DESTINATION share/pairinteraction/calc)
install(TARGETS pairinteraction-complex RUNTIME DESTINATION share/pairinteraction/calc)
install(TARGETS matrixelements RUNTIME DESTINATION share/pairinteraction/calc)

