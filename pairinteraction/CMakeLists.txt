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

target_compile_features(pireal PRIVATE cxx_std_17)
set_target_properties(pireal PROPERTIES CXX_EXTENSIONS OFF)

target_compile_definitions(pireal PUBLIC $<$<CXX_COMPILER_ID:MSVC>:NOMINMAX WIN32_LEAN_AND_MEAN _USE_MATH_DEFINES>)

add_executable(pairinteraction-real    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

target_compile_features(pairinteraction-real PRIVATE cxx_std_17)
set_target_properties(pairinteraction-real PROPERTIES CXX_EXTENSIONS OFF)

target_link_libraries(pairinteraction-real    pireal)

# Add current directory to search path

target_include_directories(pireal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(pireal PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# Set up clang-tidy

if(WITH_CLANG_TIDY)
  set(CLANG_TIDY_COMMAND_LINE "${CLANG_TIDY_EXE}" "-config=" "-warnings-as-errors=*")
  set_target_properties(pireal PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND_LINE}")
endif()

# Generate coverage report

if( WITH_COVERAGE )
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the build type" FORCE)
  target_compile_options(pireal PUBLIC $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:--coverage -fprofile-arcs -ftest-coverage>)
  target_link_libraries(pireal PUBLIC gcov)
endif( )

# Threads

find_package(Threads REQUIRED)
target_link_libraries(pireal PUBLIC Threads::Threads)

# Use OpenMP

find_package(OpenMP COMPONENTS CXX)

if (OpenMP_CXX_FOUND AND NOT WITH_CLANG_TIDY)
  target_link_libraries(pireal PUBLIC OpenMP::OpenMP_CXX)
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
set_target_properties(wignerSymbols PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(pireal PRIVATE wignerSymbols)

# Eigen
find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(pireal PUBLIC Eigen3::Eigen)

# Find MKL
set(CMAKE_C_COMPILER_ID "${CMAKE_CXX_COMPILER_ID}") # Needed by MKLConfig to find gnu_thread
set(MKL_ARCH "intel64" CACHE STRING "MKL architecture")
set(MKL_LINK "dynamic" CACHE STRING "MKL linkage")
set(MKL_THREADING "sequential" CACHE STRING "MKL threading backend")
set(MKL_INTERFACE "lp64" CACHE STRING "MKL interface")
find_package(MKL CONFIG)
find_package(PkgConfig)
if (NOT MKL_FOUND AND PKG_CONFIG_FOUND)
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
endif()

if(MSVC)
  # https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-windows/2023-1/contents-of-the-redist-intel64-directory.html
  file(GLOB MKL_DLLS
    # Threading layer
    "${MKL_ROOT}/redist/intel64/mkl_${MKL_THREADING}*.dll"
    # Computational layer
    "${MKL_ROOT}/redist/intel64/mkl_core*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_def*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_mc*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_avx*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_vml_def*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_vml_mc*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_vml_avx*.dll"
    "${MKL_ROOT}/redist/intel64/mkl_vml_cmpt*.dll"
    "${MKL_ROOT}/redist/intel64/libimalloc.dll")
  if(NOT MKL_DLLS)
    message(WARNING "MKL DLLs could not be located! Dynamic lookup will probably fail.")
  endif()
  if(MKL_THREADING STREQUAL "intel_thread")
    # Implementation detail of MKLConfig.cmake, the pkg-config files do not expose these variables
    if(EXISTS "${OMP_DLL_DIR}/${OMP_DLLNAME}")
      list(APPEND MKL_DLLS ${OMP_DLL_DIR}/${OMP_DLLNAME})
    else()
      message(WARNING "Intel OpenMP DLL cloud not be located! Dynamic lookup will probably fail.")
    endif()
  elseif(MKL_THREADING STREQUAL "tbb_thread")
    # FIXME: Do we have to add any additional DLLs?
  endif()
  set(MKL_DLLS_TO_COPY "${MKL_DLLS}" CACHE STRING "Intel MKL DLLs to copy into the target dir (;-list, ignored if not on Windows)")
  unset(MKL_DLLS)
endif()

function(target_copy_mkl_dlls TARGET)
  foreach(FILE ${MKL_DLLS_TO_COPY})
    add_custom_command(TARGET ${TARGET} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${FILE} $<TARGET_FILE_DIR:${TARGET}>)
  endforeach()
endfunction()

target_copy_mkl_dlls(pireal)

# Find a different lapack vendor
find_package(LAPACK)

# Determine which lapack vendor to use
if(WITH_LAPACK AND MKL_FOUND AND (NOT BLA_VENDOR OR BLA_VENDOR MATCHES "All" OR BLA_VENDOR MATCHES "^Intel"))
  # Use MKL
  if (TARGET MKL::MKL)
    if(OMP_LIBRARY)
      get_filename_component(OMP_LIBRARY_PATH "${OMP_LIBRARY}" DIRECTORY)
      get_filename_component(OMP_LIBRARY_PATH "${OMP_LIBRARY_PATH}" REALPATH)
      target_link_options(pireal PUBLIC "$<$<CXX_COMPILER_ID:Clang,GNU>:LINKER:-rpath=${OMP_LIBRARY_PATH}>")
    endif()
    target_link_libraries(pireal PUBLIC MKL::MKL)
    target_compile_definitions(pireal PUBLIC WITH_INTEL_MKL)
  else()
    target_link_libraries(pireal PUBLIC PkgConfig::MKL)
    target_compile_definitions(pireal PUBLIC WITH_INTEL_MKL)
  endif()

  # Using Intel MKL from Eigen
  message(STATUS "Using Intel MKL from Eigen")
  target_compile_definitions(pireal PUBLIC
    EIGEN_USE_MKL_ALL
    "-DMKL_Complex8=std::complex<float>"
    "-DMKL_Complex16=std::complex<double>")

elseif(WITH_LAPACK AND LAPACK_FOUND)
  # Use a different lapack vendor
  target_link_options(pireal PUBLIC ${LAPACK_LINKER_FLAGS})
  target_link_libraries(pireal PUBLIC ${LAPACK_LIBRARIES})

  # Some systems have LAPACKE in a separate library for some reason. The
  # following attempts to detect that. Static linking is preferred.
  add_library(LAPACKE_TARGET INTERFACE)
  add_library(LAPACKE::LAPACKE ALIAS LAPACKE_TARGET)
  find_library(LAPACKE_LIBRARY
    NAMES liblapacke.a liblapacke.dylib lapacke
    PATH_SUFFIXES lib
    PATHS /usr/local/opt/lapack)
  find_path(LAPACKE_INCLUDE_DIR
    NAMES lapacke.h
    PATH_SUFFIXES include
    PATHS /usr/include/openblas /usr/local/opt/lapack)

  set(CMAKE_REQUIRED_INCLUDES_SAVED ${CMAKE_REQUIRED_INCLUDES})
  set(CMAKE_REQUIRED_LIBRARIES_SAVED ${CMAKE_REQUIRED_LIBRARIES})
  if(LAPACKE_INCLUDE_DIR)
    set(CMAKE_REQUIRED_INCLUDES "${LAPACKE_INCLUDE_DIR}")
  endif()
  if(LAPACKE_LIBRARY)
    set(CMAKE_REQUIRED_LIBRARIES "${LAPACKE_LIBRARY};${LAPACK_LIBRARIES}")
  else()
    set(CMAKE_REQUIRED_LIBRARIES "${LAPACK_LIBRARIES}")
  endif()
  include(CheckSymbolExists)
  check_symbol_exists(LAPACKE_dsyev "lapacke.h" LAPACKE_WORKS)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES_SAVED})
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_SAVED})

  if(LAPACKE_WORKS)
    if(LAPACKE_LIBRARY)
      target_link_libraries(LAPACKE_TARGET INTERFACE "${LAPACKE_LIBRARY}")
    endif()
    if(LAPACKE_INCLUDE_DIR)
      target_include_directories(LAPACKE_TARGET SYSTEM INTERFACE "${LAPACKE_INCLUDE_DIR}")
    endif()
  elseif(NOT LAPACKE_WORKS)
    message(FATAL_ERROR "No working LAPACKE detected")
  endif()

  target_link_libraries(pireal PUBLIC LAPACKE::LAPACKE)

  # Using BLAS/LAPACK from Eigen
  message(STATUS "Using BLAS/LAPACK from Eigen")
  target_compile_definitions(pireal PUBLIC
    EIGEN_USE_BLAS EIGEN_USE_LAPACKE
    "-Dlapack_complex_float=std::complex<float>"
    "-Dlapack_complex_double=std::complex<double>")
endif()

find_package(fmt REQUIRED)
target_link_libraries(pireal PUBLIC $<IF:$<TARGET_EXISTS:fmt::fmt-header-only>,fmt::fmt-header-only,fmt::fmt>)

find_package(SQLite3 REQUIRED)
target_link_libraries(pireal PUBLIC SQLite::SQLite3)

find_package(Boost COMPONENTS serialization REQUIRED)
target_link_libraries(pireal PUBLIC Boost::boost Boost::serialization)
if(MSVC)
  target_link_libraries(pireal PUBLIC Boost::disable_autolinking Boost::dynamic_linking)
  target_compile_definitions(pireal PUBLIC BOOST_UUID_FORCE_AUTO_LINK)
endif()

if(NOT MSVC AND NOT APPLE)
  set(CMAKE_REQUIRED_LIBRARIES_SAVED ${CMAKE_REQUIRED_LIBRARIES})
  set(CMAKE_REQUIRED_LIBRARIES -lstdc++fs)
  check_cxx_source_compiles("int main() {}" CAN_LINK_STDCXXFS)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_SAVED})
  if(CAN_LINK_STDCXXFS)
    target_link_libraries(pireal PUBLIC stdc++fs)
  endif()
endif()

if(WITH_GSL)
  find_package(GSL REQUIRED)
  target_link_libraries(pireal PUBLIC GSL::gsl GSL::gslcblas)
  target_compile_definitions(pireal PUBLIC WITH_GSL)
endif()

# Generate interface with SWIG
if(WITH_PYTHON)
  # FIXME: The OLD behavior for policy CMP0078 will be removed from a future version
  # of CMake.
  if(POLICY CMP0078)
    cmake_policy(SET CMP0078 OLD)
  endif()

  set( DATATYPE "real" )
  configure_file( Interface.i.cmakein ${CMAKE_CURRENT_BINARY_DIR}/interface_real.i )

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
  target_copy_mkl_dlls(${SWIG_MODULE_pireal_REAL_NAME})

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
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pireal.py DESTINATION pairinteraction)
  endif()
  install(TARGETS pireal LIBRARY DESTINATION pairinteraction)
  install(TARGETS pairinteraction-real RUNTIME DESTINATION pairinteraction)

  set(bin1 \${CMAKE_INSTALL_PREFIX}/pairinteraction/libpireal.dylib)
  set(bin2 \${CMAKE_INSTALL_PREFIX}/pairinteraction/pairinteraction-real)
  if(WITH_PYTHON)
    set(bin3 \${CMAKE_INSTALL_PREFIX}/pairinteraction/_pireal.so)
  endif()

  install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} ${CMAKE_MACOSX_GOODIES_PATH}/standalone.py \${CMAKE_INSTALL_PREFIX}/pairinteraction/libraries ${bin1} ${bin2} ${bin3})")

elseif ( NOT WIN32 )

  if(WITH_PYTHON)
    set(PY_VERSION "${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}")
    set(PY_MOD_DIR "lib/python${PY_VERSION}/dist-packages/pairinteraction")

    install(TARGETS ${SWIG_MODULE_pireal_REAL_NAME}    LIBRARY DESTINATION ${PY_MOD_DIR})
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pireal.py    DESTINATION ${PY_MOD_DIR})
  endif()

  install(TARGETS pireal    LIBRARY DESTINATION lib)
  install(TARGETS pairinteraction-real    RUNTIME DESTINATION share/pairinteraction/pairinteraction)

endif( )
