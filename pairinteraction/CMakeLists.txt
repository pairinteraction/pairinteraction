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

# Set compiler options for release builds

if(CMAKE_BUILD_TYPE MATCHES RELEASE)
  target_compile_options(pireal
    PUBLIC $<$<CXX_COMPILER_ID:MSVC>:/Ox>
           $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Ofast -pipe>)
  target_compile_options(picomplex
    PUBLIC $<$<CXX_COMPILER_ID:MSVC>:/Ox>
           $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Ofast -pipe>)
endif()

# Add current directory to search path

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Generate in-memory database

file(READ databases/quantum_defects.sql QUANTUM_DEFECT_DATABASE_CONTENT)
configure_file(EmbeddedDatabase.h.in ${CMAKE_CURRENT_BINARY_DIR}/EmbeddedDatabase.h)

# Generate header containing the version informations
configure_file(version.h.in ${CMAKE_CURRENT_BINARY_DIR}/version.h)

# Add the executables

file(GLOB pairinteraction_SRCS *.h *.cpp)
list(REMOVE_ITEM pairinteraction_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

add_library(pireal SHARED ${pairinteraction_SRCS})
add_library(picomplex SHARED ${pairinteraction_SRCS})
target_compile_definitions(picomplex PUBLIC USE_COMPLEX)

add_executable(pairinteraction-real    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable(pairinteraction-complex ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

target_link_libraries(pairinteraction-real    pireal)
target_link_libraries(pairinteraction-complex picomplex)

if(MSVC)
  add_definitions(/W3 /DNOMINMAX /DWIN32_LEAN_AND_MEAN /D_USE_MATH_DEFINES)
  #set_source_files_properties(${pairinteraction_SRCS} PROPERTIES COMPILE_FLAGS "/W3 /DNOMINMAX /DWIN32_LEAN_AND_MEAN")
else()
  set_source_files_properties(${pairinteraction_SRCS} PROPERTIES COMPILE_FLAGS "-Wall -Wextra -pedantic -Werror")
endif()

# Set up clang-tidy

if(WITH_CLANG_TIDY)
  set(CLANG_TIDY_COMMAND_LINE "${CLANG_TIDY_EXE}" "-config=" "-warnings-as-errors=*")
  set_target_properties(pireal PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND_LINE}")
  set_target_properties(picomplex PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND_LINE}")
endif()

# Generate coverage report

if( WITH_COVERAGE )
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the build type" FORCE)
  target_compile_options(pireal PUBLIC $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-g -O0 --coverage -fprofile-arcs -ftest-coverage>)
  target_compile_options(picomplex PUBLIC $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-g -O0 --coverage -fprofile-arcs -ftest-coverage>)
  target_link_libraries(pireal PUBLIC gcov)
  target_link_libraries(picomplex PUBLIC gcov)
endif( )

# Threads

find_package(Threads REQUIRED)
target_link_libraries(pireal PUBLIC Threads::Threads)
target_link_libraries(picomplex PUBLIC Threads::Threads)

# Use OpenMP

find_package(OpenMP COMPONENTS CXX)

if (OpenMP_CXX_FOUND AND NOT WITH_CLANG_TIDY)
  target_link_libraries(pireal PUBLIC OpenMP::OpenMP_CXX)
  target_link_libraries(picomplex PUBLIC OpenMP::OpenMP_CXX)
else()
  message(WARNING "Because OpenMP was not found, calculations are not be parallelized.")
endif()

# Add external libraries

# Wigner Symbols
include(FetchContent)
FetchContent_Declare(
  wignerSymbols
  GIT_REPOSITORY https://github.com/joeydumont/wignerSymbols.git
  GIT_TAG        337dd7e634f5990f432bb53c7a69ad62c009000a
)
FetchContent_GetProperties(wignerSymbols)
if(NOT wignersymbols_POPULATED)
  FetchContent_Populate(wignerSymbols)
endif()
add_library(wignerSymbols STATIC ${wignersymbols_SOURCE_DIR}/src/wignerSymbols-cpp.cpp)
target_include_directories(pireal SYSTEM PUBLIC ${wignersymbols_SOURCE_DIR}/include)
target_include_directories(picomplex SYSTEM PUBLIC ${wignersymbols_SOURCE_DIR}/include)
set_target_properties(wignerSymbols PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(pireal PRIVATE wignerSymbols)
target_link_libraries(picomplex PRIVATE wignerSymbols)

# Include Eigen
find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(pireal PUBLIC Eigen3::Eigen)
target_link_libraries(picomplex PUBLIC Eigen3::Eigen)

find_package(SQLite3 REQUIRED)
target_link_libraries(pireal PUBLIC SQLite::SQLite3)
target_link_libraries(picomplex PUBLIC SQLite::SQLite3)

find_package(Boost COMPONENTS filesystem system serialization REQUIRED)
target_link_libraries(pireal PUBLIC Boost::boost Boost::filesystem Boost::system Boost::serialization)
target_link_libraries(picomplex PUBLIC Boost::boost Boost::filesystem Boost::system Boost::serialization)
if(MSVC)
  target_link_libraries(pireal PUBLIC Boost::disable_autolinking Boost::dynamic_linking)
  target_link_libraries(picomplex PUBLIC Boost::disable_autolinking Boost::dynamic_linking)
  target_compile_definitions(pireal PUBLIC BOOST_UUID_FORCE_AUTO_LINK)
  target_compile_definitions(picomplex PUBLIC BOOST_UUID_FORCE_AUTO_LINK)
endif()

if(WITH_GSL)
  find_package(GSL REQUIRED)
  target_link_libraries(pireal PUBLIC GSL::gsl GSL::gslcblas)
  target_link_libraries(picomplex PUBLIC GSL::gsl GSL::gslcblas)
  target_compile_definitions(pireal PUBLIC WITH_GSL)
  target_compile_definitions(picomplex PUBLIC WITH_GSL)
endif()

find_package(MKL)
if (MKL_FOUND)
  target_include_directories(pireal SYSTEM PUBLIC ${MKL_INCLUDE_DIR})
  target_include_directories(picomplex SYSTEM PUBLIC ${MKL_INCLUDE_DIR})
  target_link_libraries(pireal PUBLIC ${MKL_LIBRARY})
  target_link_libraries(picomplex PUBLIC ${MKL_LIBRARY})
  target_compile_definitions(pireal PUBLIC WITH_INTEL_MKL)
  target_compile_definitions(picomplex PUBLIC WITH_INTEL_MKL)
endif()

# Generate interface with SWIG
if(WITH_PYTHON)
  if(POLICY CMP0078)
    cmake_policy(SET CMP0078 OLD)
  endif()

  set( DATATYPE "real" )
  configure_file( Interface.i.cmakein ${CMAKE_CURRENT_BINARY_DIR}/interface_real.i )
  set( DATATYPE "complex" )
  configure_file( Interface.i.cmakein ${CMAKE_CURRENT_BINARY_DIR}/interface_complex.i )
  
  file( WRITE ${CMAKE_CURRENT_BINARY_DIR}/__init__.py "" )

  find_package(SWIG 3.0 COMPONENTS python REQUIRED)
  if(${CMAKE_VERSION} VERSION_LESS "3.18.0") 
    find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
  else()
    find_package(Python3 COMPONENTS Interpreter Development.Module NumPy REQUIRED)
  endif()

  if(MSVC)
    string(REPLACE "python3${Python3_VERSION_MINOR}.lib" "python3.lib" Python3_LIBRARIES_GENERALIZED ${Python3_LIBRARIES})
    message(STATUS "Unpatched Python3 libraries: ${Python3_LIBRARIES}")
    message(STATUS "  Patched Python3 libraries: ${Python3_LIBRARIES_GENERALIZED}")
  endif()

  include(${SWIG_USE_FILE})
  set(UseSWIG_TARGET_NAME_PREFERENCE LEGACY)

  set_property(SOURCE ${CMAKE_CURRENT_BINARY_DIR}/interface_real.i PROPERTY CPLUSPLUS ON)
  swig_add_library(pireal LANGUAGE python SOURCES ${CMAKE_CURRENT_BINARY_DIR}/interface_real.i)
  set_target_properties(${SWIG_MODULE_pireal_REAL_NAME} PROPERTIES
    SWIG_COMPILE_OPTIONS -py3
    SWIG_GENERATED_INCLUDE_DIRECTORIES "${Python3_INCLUDE_DIRS};${Python3_NumPy_INCLUDE_DIRS}"
    SWIG_GENERATED_COMPILE_OPTIONS $<$<CXX_COMPILER_ID:MSVC>:/bigobj>)
  target_link_options(${SWIG_MODULE_pireal_REAL_NAME} PUBLIC
    $<$<CXX_COMPILER_ID:AppleClang>:-undefined dynamic_lookup>
    $<$<CXX_COMPILER_ID:Clang,GNU>:-Wl,--unresolved-symbols=ignore-all>
    $<$<CXX_COMPILER_ID:MSVC>:/NODEFAULTLIB:python3${Python3_VERSION_MINOR}.lib /DEFAULTLIB:python3.lib /FORCE:UNRESOLVED>)
  swig_link_libraries(pireal PUBLIC pireal $<$<CXX_COMPILER_ID:MSVC>:${Python3_LIBRARIES_GENERALIZED}>)

  set_property(SOURCE ${CMAKE_CURRENT_BINARY_DIR}/interface_complex.i PROPERTY CPLUSPLUS ON)
  swig_add_library(picomplex LANGUAGE python SOURCES ${CMAKE_CURRENT_BINARY_DIR}/interface_complex.i)
  set_target_properties(${SWIG_MODULE_picomplex_REAL_NAME} PROPERTIES
    SWIG_COMPILE_OPTIONS "-py3;-DUSE_COMPLEX"
    SWIG_GENERATED_COMPILE_DEFINITIONS USE_COMPLEX
    SWIG_GENERATED_INCLUDE_DIRECTORIES "${Python3_INCLUDE_DIRS};${Python3_NumPy_INCLUDE_DIRS}"
    SWIG_GENERATED_COMPILE_OPTIONS $<$<CXX_COMPILER_ID:MSVC>:/bigobj>)
  target_link_options(${SWIG_MODULE_picomplex_REAL_NAME} PUBLIC
    $<$<CXX_COMPILER_ID:AppleClang>:-undefined dynamic_lookup>
    $<$<CXX_COMPILER_ID:Clang,GNU>:-Wl,--unresolved-symbols=ignore-all>
    $<$<CXX_COMPILER_ID:MSVC>:/NODEFAULTLIB:python3${Python3_VERSION_MINOR}.lib /DEFAULTLIB:python3.lib /FORCE:UNRESOLVED>)
  swig_link_libraries(picomplex PUBLIC picomplex $<$<CXX_COMPILER_ID:MSVC>:${Python3_LIBRARIES_GENERALIZED}>)

endif()

# Include test cases

add_subdirectory( unit_test )

# CPack

if( WITH_DMG )
  find_package(Python3 COMPONENTS Interpreter REQUIRED)

  if(WITH_PYTHON)
    install(TARGETS ${SWIG_MODULE_pireal_REAL_NAME} LIBRARY DESTINATION ${LIBNAME})
    install(TARGETS ${SWIG_MODULE_picomplex_REAL_NAME} LIBRARY DESTINATION ${LIBNAME})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pireal.py DESTINATION ${LIBNAME})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/picomplex.py DESTINATION ${LIBNAME})
  endif()
  install(TARGETS pireal LIBRARY DESTINATION ${LIBNAME})
  install(TARGETS picomplex LIBRARY DESTINATION ${LIBNAME})
  install(TARGETS pairinteraction-real RUNTIME DESTINATION ${LIBNAME})
  install(TARGETS pairinteraction-complex RUNTIME DESTINATION ${LIBNAME})

  set(bin1 \${CMAKE_INSTALL_PREFIX}/${LIBNAME}/libpireal.dylib)
  set(bin2 \${CMAKE_INSTALL_PREFIX}/${LIBNAME}/libpicomplex.dylib)
  set(bin3 \${CMAKE_INSTALL_PREFIX}/${LIBNAME}/pairinteraction-real)
  set(bin4 \${CMAKE_INSTALL_PREFIX}/${LIBNAME}/pairinteraction-complex)
  if(WITH_PYTHON)
    set(bin5 \${CMAKE_INSTALL_PREFIX}/${LIBNAME}/_pireal.so)
    set(bin6 \${CMAKE_INSTALL_PREFIX}/${LIBNAME}/_picomplex.so)
  endif()
  
  install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} ${CMAKE_MACOSX_GOODIES_PATH}/standalone.py \${CMAKE_INSTALL_PREFIX}/${LIBNAME}/libraries ${bin1} ${bin2} ${bin3} ${bin4} ${bin5} ${bin6})")

elseif ( NOT WIN32 )

  if(WITH_PYTHON)
    set(PY_VERSION "${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}")
    set(PY_MOD_DIR "lib/python${PY_VERSION}/dist-packages/${LIBNAME}")

    install(TARGETS ${SWIG_MODULE_pireal_REAL_NAME}    LIBRARY DESTINATION ${PY_MOD_DIR})
    install(TARGETS ${SWIG_MODULE_picomplex_REAL_NAME} LIBRARY DESTINATION ${PY_MOD_DIR})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pireal.py    DESTINATION ${PY_MOD_DIR})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/picomplex.py DESTINATION ${PY_MOD_DIR})
  endif()

  install(TARGETS pireal    LIBRARY DESTINATION lib)
  install(TARGETS picomplex LIBRARY DESTINATION lib)
  install(TARGETS pairinteraction-real    RUNTIME DESTINATION share/pairinteraction/${LIBNAME})
  install(TARGETS pairinteraction-complex RUNTIME DESTINATION share/pairinteraction/${LIBNAME})

endif( )
