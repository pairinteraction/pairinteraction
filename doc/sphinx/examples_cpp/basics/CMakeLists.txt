# Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
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

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(basics CXX)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Ofast -pipe")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add executable

add_executable(basics basics.cpp)

# Enable OpenMP if supported

find_package(OpenMP COMPONENTS CXX)
if (OpenMP_CXX_FOUND})
    target_link_libraries(basics PUBLIC OpenMP::OpenMP_CXX)
endif()

# Link the pairinteraction library and its dependencies # TODO create a FindPairinteraction.cmake file

include_directories(SYSTEM TODO/pairinteraction/eigen)

#target_compile_definitions(basics PRIVATE USE_COMPLEX) # if used, libpicomplex.dylib has to be linked
include_directories(TODO/pairinteraction)
list(APPEND LIBRARIES "TODO/pairinteraction/build/pairinteraction/libpireal.dylib")

# Link boost filesystem

find_package(Boost COMPONENTS filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
list(APPEND LIBRARIES ${Boost_LIBRARIES})

# Link libraries

target_link_libraries(basics ${LIBRARIES})
