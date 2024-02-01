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
list(REMOVE_ITEM pairinteraction_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/binding.cpp)
list(REMOVE_ITEM pairinteraction_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

add_library(pairinteraction SHARED ${pairinteraction_SRCS})
set_target_properties(pairinteraction PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_compile_features(pairinteraction PRIVATE cxx_std_17)
set_target_properties(pairinteraction PROPERTIES CXX_EXTENSIONS OFF)

target_compile_definitions(pairinteraction PUBLIC $<$<CXX_COMPILER_ID:MSVC>:NOMINMAX WIN32_LEAN_AND_MEAN _USE_MATH_DEFINES>)

add_executable(pairinteraction-backend-deprecated ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

target_compile_features(pairinteraction-backend-deprecated PRIVATE cxx_std_17)
set_target_properties(pairinteraction-backend-deprecated PROPERTIES
  CXX_EXTENSIONS OFF
  INSTALL_RPATH $<$<BOOL:${SKBUILD}>:$<$<PLATFORM_ID:Linux>:$ORIGIN>$<$<PLATFORM_ID:Darwin>:@loader_path>>)

target_link_libraries(pairinteraction-backend-deprecated pairinteraction)

# Add current directory to search path

target_include_directories(pairinteraction PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(pairinteraction PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# Generate coverage report

if( WITH_COVERAGE )
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the build type" FORCE)
  target_compile_options(pairinteraction PUBLIC $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:--coverage -fprofile-arcs -ftest-coverage>)
  target_link_libraries(pairinteraction PUBLIC gcov)
endif( )

# Threads

find_package(Threads REQUIRED)
target_link_libraries(pairinteraction PUBLIC Threads::Threads)

# Use OpenMP

find_package(OpenMP COMPONENTS CXX)

if (OpenMP_CXX_FOUND AND NOT WITH_CLANG_TIDY)
  target_link_libraries(pairinteraction PUBLIC OpenMP::OpenMP_CXX)
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
target_include_directories(pairinteraction SYSTEM PUBLIC ${wignersymbols_SOURCE_DIR}/include)
set_target_properties(wignerSymbols PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CXX_CLANG_TIDY "")
target_link_libraries(pairinteraction PRIVATE wignerSymbols)

# Eigen
find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(pairinteraction PUBLIC Eigen3::Eigen)

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

target_copy_mkl_dlls(pairinteraction)

# Find a different lapack vendor
find_package(LAPACK)

# Determine which lapack vendor to use
if(WITH_LAPACK AND MKL_FOUND AND (NOT BLA_VENDOR OR BLA_VENDOR MATCHES "All" OR BLA_VENDOR MATCHES "^Intel"))
  # Use MKL
  if (TARGET MKL::MKL)
    if(OMP_LIBRARY)
      get_filename_component(OMP_LIBRARY_PATH "${OMP_LIBRARY}" DIRECTORY)
      get_filename_component(OMP_LIBRARY_PATH "${OMP_LIBRARY_PATH}" REALPATH)
      target_link_options(pairinteraction PUBLIC "$<$<CXX_COMPILER_ID:Clang,GNU>:LINKER:-rpath=${OMP_LIBRARY_PATH}>")
    endif()
  else()
    add_library(MKL::MKL ALIAS PkgConfig::MKL)
  endif()
  target_link_libraries(pairinteraction PUBLIC MKL::MKL)
  target_compile_definitions(pairinteraction PUBLIC WITH_INTEL_MKL)

  # Using Intel MKL from Eigen
  message(STATUS "Using Intel MKL from Eigen")
  target_compile_definitions(pairinteraction PUBLIC
    EIGEN_USE_MKL_ALL
    "-DMKL_Complex8=std::complex<float>"
    "-DMKL_Complex16=std::complex<double>")

elseif(WITH_LAPACK AND LAPACK_FOUND)
  # Use a different lapack vendor
  if(NOT TARGET LAPACK::LAPACK)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    set_target_properties(LAPACK::LAPACK PROPERTIES
      INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
    set_target_properties(LAPACK::LAPACK PROPERTIES
      INTERFACE_LINK_OPTIONS "${LAPACK_LINKER_FLAGS}")
  endif()
  target_link_libraries(pairinteraction PUBLIC LAPACK::LAPACK)

  # Some systems have LAPACKE in a separate library for some reason. The
  # following attempts to detect that. Static linking is preferred.
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
    add_library(LAPACKE::LAPACKE INTERFACE IMPORTED)
    if(LAPACKE_LIBRARY)
      set_target_properties(LAPACKE::LAPACKE PROPERTIES
        INTERFACE_LINK_LIBRARIES "${LAPACKE_LIBRARY}")
    endif()
    if(LAPACKE_INCLUDE_DIR)
      set_target_properties(LAPACKE::LAPACKE PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${LAPACKE_INCLUDE_DIR}")
    endif()
  elseif(NOT LAPACKE_WORKS)
    message(FATAL_ERROR "No working LAPACKE detected")
  endif()

  target_link_libraries(pairinteraction PUBLIC LAPACKE::LAPACKE)

  # Using BLAS/LAPACK from Eigen
  message(STATUS "Using BLAS/LAPACK from Eigen")
  target_compile_definitions(pairinteraction PUBLIC
    EIGEN_USE_BLAS EIGEN_USE_LAPACKE
    "-Dlapack_complex_float=std::complex<float>"
    "-Dlapack_complex_double=std::complex<double>")
endif()

find_package(fmt REQUIRED)
target_link_libraries(pairinteraction PUBLIC $<IF:$<TARGET_EXISTS:fmt::fmt-header-only>,fmt::fmt-header-only,fmt::fmt>)

# We do not directly use "find_package(cereal CONFIG)" because of missing CMake support files in Ubuntu
# https://bugs.launchpad.net/ubuntu/+source/libcereal/+bug/1978462
add_library(CEREAL_TARGET INTERFACE)
add_library(cereal::cereal ALIAS CEREAL_TARGET)
find_path(CEREAL_INCLUDE_DIR cereal/cereal.hpp)
if(NOT CEREAL_INCLUDE_DIR)
  message(FATAL_ERROR "Cereal could not be found")
endif()
target_include_directories(CEREAL_TARGET SYSTEM INTERFACE "${CEREAL_INCLUDE_DIR}")
target_link_libraries(pairinteraction PUBLIC cereal::cereal)

find_package(SQLite3 REQUIRED)
target_link_libraries(pairinteraction PUBLIC SQLite::SQLite3)

find_package(Boost REQUIRED)
target_link_libraries(pairinteraction PUBLIC Boost::boost)
if(MSVC)
  target_link_libraries(pairinteraction PUBLIC Boost::disable_autolinking Boost::dynamic_linking)
  target_compile_definitions(pairinteraction PUBLIC BOOST_UUID_FORCE_AUTO_LINK)
endif()

if(NOT MSVC AND NOT APPLE)
  set(CMAKE_REQUIRED_LIBRARIES_SAVED ${CMAKE_REQUIRED_LIBRARIES})
  set(CMAKE_REQUIRED_LIBRARIES -lstdc++fs)
  check_cxx_source_compiles("int main() {}" CAN_LINK_STDCXXFS)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_SAVED})
  if(CAN_LINK_STDCXXFS)
    target_link_libraries(pairinteraction PUBLIC stdc++fs)
  endif()
endif()

if(WITH_GSL)
  find_package(GSL REQUIRED)
  target_link_libraries(pairinteraction PUBLIC GSL::gsl GSL::gslcblas)
  target_compile_definitions(pairinteraction PUBLIC WITH_GSL)
endif()

# Generate interface with SWIG
if(WITH_PYTHON)
  file( WRITE ${CMAKE_CURRENT_BINARY_DIR}/__init__.py "from .binding import *" )
  file( WRITE ${CMAKE_CURRENT_BINARY_DIR}/pireal.py "from .binding import *\nSystemOne,SystemTwo=SystemOneReal,SystemTwoReal" )
  file( WRITE ${CMAKE_CURRENT_BINARY_DIR}/picomplex.py "from .binding import *\nSystemOne,SystemTwo=SystemOneComplex,SystemTwoComplex" )

  if(${CMAKE_VERSION} VERSION_LESS "3.18.0")
    find_package(Python REQUIRED
      COMPONENTS Interpreter Development NumPy
      OPTIONAL_COMPONENTS Development.SABIModule)
  else()
    find_package(Python REQUIRED
      COMPONENTS Interpreter Development.Module NumPy
      OPTIONAL_COMPONENTS Development.SABIModule)
  endif()

  FetchContent_Declare(nanobind
    GIT_REPOSITORY https://github.com/wjakob/nanobind.git
    GIT_TAG        72831bb)
  FetchContent_GetProperties(nanobind)
  if(NOT nanobind_POPULATED)
    FetchContent_Populate(nanobind)
  endif()
  list(APPEND CMAKE_PREFIX_PATH ${nanobind_SOURCE_DIR})
  find_package(nanobind CONFIG REQUIRED)

  nanobind_add_module(binding STABLE_ABI NB_STATIC binding.cpp)
  set_target_properties(nanobind-static PROPERTIES CXX_CLANG_TIDY "")

  target_compile_features(binding PRIVATE cxx_std_17)
  set_target_properties(binding PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_CLANG_TIDY "")
  target_link_libraries(binding PUBLIC pairinteraction Python::NumPy)
  target_copy_mkl_dlls(binding)
  get_property(SWIG_MODULE_binding_SUPPORT_FILES TARGET binding PROPERTY SWIG_SUPPORT_FILES)

endif()

# CPack

if( WITH_DMG )
  find_package(Python COMPONENTS Interpreter REQUIRED)

  if(WITH_PYTHON)
    install(TARGETS binding LIBRARY DESTINATION pairinteraction)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pireal.py DESTINATION pairinteraction)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/picomplex.py DESTINATION pairinteraction)
  endif()
  install(TARGETS pairinteraction LIBRARY DESTINATION pairinteraction)
  install(TARGETS pairinteraction-backend-deprecated RUNTIME DESTINATION pairinteraction)

  set(bin1 \${CMAKE_INSTALL_PREFIX}/pairinteraction/libpairinteraction.dylib)
  set(bin2 \${CMAKE_INSTALL_PREFIX}/pairinteraction/pairinteraction-backend-deprecated)
  if(WITH_PYTHON)
    set(bin3 \${CMAKE_INSTALL_PREFIX}/pairinteraction/_binding.so)
  endif()

  install(CODE "execute_process(COMMAND ${Python_EXECUTABLE} ${CMAKE_MACOSX_GOODIES_PATH}/standalone.py \${CMAKE_INSTALL_PREFIX}/pairinteraction/libraries ${bin1} ${bin2} ${bin3})")

elseif(SKBUILD)
  install(TARGETS binding LIBRARY DESTINATION ${SKBUILD_PLATLIB_DIR}/pairinteraction)
  install(TARGETS pairinteraction LIBRARY DESTINATION ${SKBUILD_PLATLIB_DIR}/pairinteraction)
  install(TARGETS pairinteraction-backend-deprecated RUNTIME DESTINATION ${SKBUILD_PLATLIB_DIR}/pairinteraction)
  install(FILES ${SWIG_MODULE_binding_SUPPORT_FILES} DESTINATION ${SKBUILD_PLATLIB_DIR}/pairinteraction)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pireal.py DESTINATION ${SKBUILD_PLATLIB_DIR}/pairinteraction)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/picomplex.py DESTINATION ${SKBUILD_PLATLIB_DIR}/pairinteraction)
  install(FILES ${MKL_DLLS_TO_COPY} DESTINATION ${SKBUILD_PLATLIB_DIR}/pairinteraction)
elseif(WITH_PYTHON)
  # Hopefully better in the future... https://gitlab.kitware.com/cmake/cmake/-/issues/24213
  set(PY_MOD_DIR "lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/pairinteraction")
  install(TARGETS binding LIBRARY DESTINATION ${PY_MOD_DIR})
  install(TARGETS pairinteraction LIBRARY DESTINATION lib)
  install(TARGETS pairinteraction-backend-deprecated RUNTIME DESTINATION share/pairinteraction/pairinteraction)
  install(FILES ${SWIG_MODULE_binding_SUPPORT_FILES} DESTINATION ${PY_MOD_DIR})
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pireal.py DESTINATION ${PY_MOD_DIR})
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/picomplex.py DESTINATION ${PY_MOD_DIR})
else()
  install(TARGETS pairinteraction LIBRARY DESTINATION lib)
  install(TARGETS pairinteraction-backend-deprecated RUNTIME DESTINATION share/pairinteraction/pairinteraction)
endif()
