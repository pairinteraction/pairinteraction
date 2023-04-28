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

# Generate in-memory database

file(READ databases/quantum_defects.sql QUANTUM_DEFECT_DATABASE_CONTENT)
configure_file(EmbeddedDatabase.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/EmbeddedDatabase.hpp)

# Generate header containing the version informations
configure_file(version.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/version.hpp)

# Add the executables

file(GLOB pairinteraction_SRCS *.h *.cpp)
list(REMOVE_ITEM pairinteraction_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

add_library(pireal SHARED ${pairinteraction_SRCS})
add_library(picomplex SHARED ${pairinteraction_SRCS})
target_compile_definitions(picomplex PUBLIC USE_COMPLEX)

target_compile_features(pireal PRIVATE cxx_std_17)
set_target_properties(pireal PROPERTIES CXX_EXTENSIONS OFF)
target_compile_features(picomplex PRIVATE cxx_std_17)
set_target_properties(picomplex PROPERTIES CXX_EXTENSIONS OFF)

target_compile_definitions(pireal PUBLIC $<$<CXX_COMPILER_ID:MSVC>:NOMINMAX WIN32_LEAN_AND_MEAN _USE_MATH_DEFINES>)
target_compile_definitions(picomplex PUBLIC $<$<CXX_COMPILER_ID:MSVC>:NOMINMAX WIN32_LEAN_AND_MEAN _USE_MATH_DEFINES>)

add_executable(pairinteraction-real    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable(pairinteraction-complex ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

target_compile_features(pairinteraction-real PRIVATE cxx_std_17)
set_target_properties(pairinteraction-real PROPERTIES CXX_EXTENSIONS OFF)
target_compile_features(pairinteraction-complex PRIVATE cxx_std_17)
set_target_properties(pairinteraction-complex PROPERTIES CXX_EXTENSIONS OFF)

target_link_libraries(pairinteraction-real    pireal)
target_link_libraries(pairinteraction-complex picomplex)

# Add current directory to search path

target_include_directories(pireal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(picomplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(pireal PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(picomplex PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# Set up clang-tidy

if(WITH_CLANG_TIDY)
  set(CLANG_TIDY_COMMAND_LINE "${CLANG_TIDY_EXE}" "-config=" "-warnings-as-errors=*")
  set_target_properties(pireal PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND_LINE}")
  set_target_properties(picomplex PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND_LINE}")
endif()

# Generate coverage report

if( WITH_COVERAGE )
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the build type" FORCE)
  target_compile_options(pireal PUBLIC $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:--coverage -fprofile-arcs -ftest-coverage>)
  target_compile_options(picomplex PUBLIC $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:--coverage -fprofile-arcs -ftest-coverage>)
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

set(MKL_ARCH "intel64" CACHE STRING "MKL architecture")
set(MKL_LINK "dynamic" CACHE STRING "MKL linkage")
set(MKL_THREADING "sequential" CACHE STRING "MKL threading backend")
set(MKL_INTERFACE "lp64" CACHE STRING "MKL interface")
find_package(MKL CONFIG)
find_package(PkgConfig)
if(MKL_FOUND)
  target_link_libraries(pireal PUBLIC MKL::MKL)
  target_link_libraries(picomplex PUBLIC MKL::MKL)
  target_compile_definitions(pireal PUBLIC WITH_INTEL_MKL)
  target_compile_definitions(picomplex PUBLIC WITH_INTEL_MKL)
elseif(PKG_CONFIG_FOUND)
  if(MKL_LINK STREQUAL "sdl")
    pkg_check_modules(MKL mkl-sdl IMPORTED_TARGET)
  elseif(MKL_THREADING STREQUAL "sequential")
    pkg_check_modules(MKL mkl-${MKL_LINK}-${MKL_INTERFACE}-seq IMPORTED_TARGET)
  elseif(MKL_THREADING STREQUAL "intel_thread")
    pkg_check_modules(MKL mkl-${MKL_LINK}-${MKL_INTERFACE}-iomp IMPORTED_TARGET)
  elseif(MKL_THREADING STREQUAL "gnu_thread")
    pkg_check_modules(MKL mkl-${MKL_LINK}-${MKL_INTERFACE}-gomp IMPORTED_TARGET)
  elseif(MKL_THREADING STREQUAL "tbb_thread")
    pkg_check_modules(MKL mkl-${MKL_LINK}-${MKL_INTERFACE}-tbb IMPORTED_TARGET)
  endif()
  if(MKL_FOUND)
    target_link_libraries(pireal PUBLIC PkgConfig::MKL)
    target_link_libraries(picomplex PUBLIC PkgConfig::MKL)
    target_compile_definitions(pireal PUBLIC WITH_INTEL_MKL)
    target_compile_definitions(picomplex PUBLIC WITH_INTEL_MKL)
  endif()
endif()

find_package(BLAS)
find_package(LAPACK)
if(WITH_LAPACK AND LAPACK_FOUND AND NOT MKL_FOUND)
  target_link_options(pireal PUBLIC ${LAPACK_LINKER_FLAGS})
  target_link_options(picomplex PUBLIC ${LAPACK_LINKER_FLAGS})
  target_link_libraries(pireal PUBLIC ${LAPACK_LIBRARIES})
  target_link_libraries(picomplex PUBLIC ${LAPACK_LIBRARIES})

  # TODO: Ubuntu and macOS have this in a separate library for some reason.
  # FIXME: Needs a more robust way to detection and linking.
  find_library(LAPACKE_LIBRARY NAMES liblapacke.a)

  set(CMAKE_REQUIRED_LIBRARIES_SAVED ${CMAKE_REQUIRED_LIBRARIES})
  if(LAPACKE_LIBRARY)
    set(CMAKE_REQUIRED_LIBRARIES "${LAPACKE_LIBRARY};${LAPACK_LIBRARIES}")
  else()
    set(CMAKE_REQUIRED_LIBRARIES "${LAPACK_LIBRARIES}")
  endif()
  check_function_exists(LAPACKE_dsyev LAPACKE_WORKS)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_SAVED})

  if(LAPACKE_WORKS AND LAPACKE_LIBRARY)
    target_link_libraries(pireal PUBLIC ${LAPACKE_LIBRARY})
    target_link_libraries(picomplex PUBLIC ${LAPACKE_LIBRARY})
  elseif(NOT LAPACKE_WORKS)
    message(FATAL_ERROR "No working LAPACKE detected")
  endif()
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(pireal PUBLIC Eigen3::Eigen)
target_link_libraries(picomplex PUBLIC Eigen3::Eigen)
if(WITH_LAPACK AND MKL_FOUND)
  target_compile_definitions(pireal PUBLIC
    EIGEN_USE_MKL_ALL
    "-DMKL_Complex8=std::complex<float>"
    "-DMKL_Complex16=std::complex<double>")
  target_compile_definitions(picomplex PUBLIC
    EIGEN_USE_MKL_ALL
    "-DMKL_Complex8=std::complex<float>"
    "-DMKL_Complex16=std::complex<double>")
elseif(WITH_LAPACK AND LAPACK_FOUND)
  target_compile_definitions(pireal PUBLIC
    EIGEN_USE_BLAS EIGEN_USE_LAPACKE
    "-Dlapack_complex_float=std::complex<float>"
    "-Dlapack_complex_double=std::complex<double>")
  target_compile_definitions(picomplex PUBLIC
    EIGEN_USE_BLAS EIGEN_USE_LAPACKE
    "-Dlapack_complex_float=std::complex<float>"
    "-Dlapack_complex_double=std::complex<double>")
endif()

find_package(fmt REQUIRED)
target_link_libraries(pireal PUBLIC $<IF:$<TARGET_EXISTS:fmt::fmt-header-only>,fmt::fmt-header-only,fmt::fmt>)
target_link_libraries(picomplex PUBLIC $<IF:$<TARGET_EXISTS:fmt::fmt-header-only>,fmt::fmt-header-only,fmt::fmt>)

find_package(SQLite3 REQUIRED)
target_link_libraries(pireal PUBLIC SQLite::SQLite3)
target_link_libraries(picomplex PUBLIC SQLite::SQLite3)

find_package(Boost COMPONENTS serialization REQUIRED)
target_link_libraries(pireal PUBLIC Boost::boost Boost::serialization)
target_link_libraries(picomplex PUBLIC Boost::boost Boost::serialization)
if(MSVC)
  target_link_libraries(pireal PUBLIC Boost::disable_autolinking Boost::dynamic_linking)
  target_link_libraries(picomplex PUBLIC Boost::disable_autolinking Boost::dynamic_linking)
  target_compile_definitions(pireal PUBLIC BOOST_UUID_FORCE_AUTO_LINK)
  target_compile_definitions(picomplex PUBLIC BOOST_UUID_FORCE_AUTO_LINK)
endif()

if(NOT MSVC AND NOT APPLE)
  set(CMAKE_REQUIRED_LIBRARIES_SAVED ${CMAKE_REQUIRED_LIBRARIES})
  set(CMAKE_REQUIRED_LIBRARIES -lstdc++fs)
  check_cxx_source_compiles("int main() {}" CAN_LINK_STDCXXFS)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_SAVED})
  if(CAN_LINK_STDCXXFS)
    target_link_libraries(pireal PUBLIC stdc++fs)
    target_link_libraries(picomplex PUBLIC stdc++fs)
  endif()
endif()

if(WITH_GSL)
  find_package(GSL REQUIRED)
  target_link_libraries(pireal PUBLIC GSL::gsl GSL::gslcblas)
  target_link_libraries(picomplex PUBLIC GSL::gsl GSL::gslcblas)
  target_compile_definitions(pireal PUBLIC WITH_GSL)
  target_compile_definitions(picomplex PUBLIC WITH_GSL)
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
  target_compile_features(${SWIG_MODULE_pireal_REAL_NAME} PRIVATE cxx_std_17)
  set_target_properties(${SWIG_MODULE_pireal_REAL_NAME} PROPERTIES
    CXX_EXTENSIONS OFF
    SWIG_COMPILE_OPTIONS "-Wall;-py3"
    SWIG_USE_TARGET_INCLUDE_DIRECTORIES ON
    SWIG_GENERATED_INCLUDE_DIRECTORIES "${Python3_INCLUDE_DIRS};${Python3_NumPy_INCLUDE_DIRS}"
    SWIG_GENERATED_COMPILE_OPTIONS $<$<CXX_COMPILER_ID:MSVC>:/bigobj>)
  target_link_options(${SWIG_MODULE_pireal_REAL_NAME} PUBLIC
    $<$<CXX_COMPILER_ID:AppleClang>:-undefined dynamic_lookup>
    $<$<CXX_COMPILER_ID:Clang,GNU>:-Wl,--unresolved-symbols=ignore-all>
    $<$<CXX_COMPILER_ID:MSVC>:/NODEFAULTLIB:python3${Python3_VERSION_MINOR}.lib /DEFAULTLIB:python3.lib /FORCE:UNRESOLVED>)
  swig_link_libraries(pireal PUBLIC pireal $<$<CXX_COMPILER_ID:MSVC>:${Python3_LIBRARIES_GENERALIZED}>)

  set_property(SOURCE ${CMAKE_CURRENT_BINARY_DIR}/interface_complex.i PROPERTY CPLUSPLUS ON)
  swig_add_library(picomplex LANGUAGE python SOURCES ${CMAKE_CURRENT_BINARY_DIR}/interface_complex.i)
  target_compile_features(${SWIG_MODULE_picomplex_REAL_NAME} PRIVATE cxx_std_17)
  set_target_properties(${SWIG_MODULE_picomplex_REAL_NAME} PROPERTIES
    CXX_EXTENSIONS OFF
    SWIG_COMPILE_OPTIONS "-Wall;-py3;-DUSE_COMPLEX"
    SWIG_USE_TARGET_INCLUDE_DIRECTORIES ON
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

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  add_subdirectory(unit_test)
endif()

# CPack

if( WITH_DMG )
  find_package(Python3 COMPONENTS Interpreter REQUIRED)

  if(WITH_PYTHON)
    install(TARGETS ${SWIG_MODULE_pireal_REAL_NAME} LIBRARY DESTINATION pairinteraction)
    install(TARGETS ${SWIG_MODULE_picomplex_REAL_NAME} LIBRARY DESTINATION pairinteraction)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pireal.py DESTINATION pairinteraction)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/picomplex.py DESTINATION pairinteraction)
  endif()
  install(TARGETS pireal LIBRARY DESTINATION pairinteraction)
  install(TARGETS picomplex LIBRARY DESTINATION pairinteraction)
  install(TARGETS pairinteraction-real RUNTIME DESTINATION pairinteraction)
  install(TARGETS pairinteraction-complex RUNTIME DESTINATION pairinteraction)

  set(bin1 \${CMAKE_INSTALL_PREFIX}/pairinteraction/libpireal.dylib)
  set(bin2 \${CMAKE_INSTALL_PREFIX}/pairinteraction/libpicomplex.dylib)
  set(bin3 \${CMAKE_INSTALL_PREFIX}/pairinteraction/pairinteraction-real)
  set(bin4 \${CMAKE_INSTALL_PREFIX}/pairinteraction/pairinteraction-complex)
  if(WITH_PYTHON)
    set(bin5 \${CMAKE_INSTALL_PREFIX}/pairinteraction/_pireal.so)
    set(bin6 \${CMAKE_INSTALL_PREFIX}/pairinteraction/_picomplex.so)
  endif()

  install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} ${CMAKE_MACOSX_GOODIES_PATH}/standalone.py \${CMAKE_INSTALL_PREFIX}/pairinteraction/libraries ${bin1} ${bin2} ${bin3} ${bin4} ${bin5} ${bin6})")

elseif ( NOT WIN32 )

  if(WITH_PYTHON)
    set(PY_VERSION "${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}")
    set(PY_MOD_DIR "lib/python${PY_VERSION}/dist-packages/pairinteraction")

    install(TARGETS ${SWIG_MODULE_pireal_REAL_NAME}    LIBRARY DESTINATION ${PY_MOD_DIR})
    install(TARGETS ${SWIG_MODULE_picomplex_REAL_NAME} LIBRARY DESTINATION ${PY_MOD_DIR})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pireal.py    DESTINATION ${PY_MOD_DIR})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/picomplex.py DESTINATION ${PY_MOD_DIR})
  endif()

  install(TARGETS pireal    LIBRARY DESTINATION lib)
  install(TARGETS picomplex LIBRARY DESTINATION lib)
  install(TARGETS pairinteraction-real    RUNTIME DESTINATION share/pairinteraction/pairinteraction)
  install(TARGETS pairinteraction-complex RUNTIME DESTINATION share/pairinteraction/pairinteraction)

endif( )
