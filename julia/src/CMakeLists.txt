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

find_program(Julia_EXECUTABLE NAMES julia REQUIRED)
if(NOT Julia_EXECUTABLE)
  message(FATAL_ERROR "Julia interpreter not found")
endif()
execute_process(
  COMMAND "${Julia_EXECUTABLE}" -e "using CxxWrap; print(joinpath(dirname(CxxWrap.jlcxx_path), \"cmake\", \"JlCxx\"))"
  OUTPUT_VARIABLE JlCxx_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)
find_package(JlCxx REQUIRED)

# Undo all the JlCxx shenanigans
unset(CMAKE_CXX_STANDRD)
set_property(TARGET JlCxx::cxxwrap_julia PROPERTY INTERFACE_COMPILE_OPTIONS)

add_library(pairinteraction_jl SHARED julia_hooks_pi.cpp)

target_compile_features(pairinteraction_jl PRIVATE cxx_std_17)
set_target_properties(pairinteraction_jl PROPERTIES
  CXX_EXTENSIONS OFF
  CXX_CLANG_TIDY "")

target_link_libraries(pairinteraction_jl PUBLIC pairinteraction JlCxx::cxxwrap_julia)

configure_file(PairInteraction.jl ${CMAKE_CURRENT_BINARY_DIR}/PairInteraction.jl COPYONLY)
