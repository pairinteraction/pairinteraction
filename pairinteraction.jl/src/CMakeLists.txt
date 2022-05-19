# Copyright (c) 2020 Sebastian Weber, Henri Menke, Alexander Papageorge. All rights reserved.
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

find_package(JuliaInterp REQUIRED)
execute_process(
  COMMAND "${Julia_EXECUTABLE}" -e "using CxxWrap; print(joinpath(dirname(CxxWrap.jlcxx_path), \"cmake\", \"JlCxx\"))"
  OUTPUT_VARIABLE JlCxx_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)
find_package(JlCxx REQUIRED)
set(CMAKE_CXX_STANDARD 17) # JlCxx overwrites this!

add_library(pireal_jl SHARED julia_hooks_pi.cpp)
add_library(picomplex_jl SHARED julia_hooks_pi.cpp)

target_compile_definitions(picomplex_jl PUBLIC "USE_COMPLEX")

include_directories(${CMAKE_SOURCE_DIR}/${LIBNAME})
include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/libs/eigen)
include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/libs/wignerSymbols/include)

target_link_libraries(pireal_jl pireal JlCxx::cxxwrap_julia)
target_link_libraries(picomplex_jl picomplex JlCxx::cxxwrap_julia)

configure_file(PairInteraction.jl ${CMAKE_CURRENT_BINARY_DIR}/PairInteraction.jl COPYONLY)
