# Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Find ZeroMQ
#
# ZEROMQ_LIBRARY       - the ZeroMQ library
# ZEROMQ_INCLUDE_DIR   - path including ZeroMQ.h

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include (FindPackageHandleStandardArgs)

# Get include dir
find_path( ZEROMQ_INCLUDE_DIR NAMES zmq.h )

# Get version, see https://github.com/zeromq/azmq/blob/master/config/FindZeroMQ.cmake
set(_ZeroMQ_H ${ZEROMQ_INCLUDE_DIR}/zmq.h)

function(_zmqver_EXTRACT _ZeroMQ_VER_COMPONENT _ZeroMQ_VER_OUTPUT)
    set(CMAKE_MATCH_1 "0")
    set(_ZeroMQ_expr "^[ \\t]*#define[ \\t]+${_ZeroMQ_VER_COMPONENT}[ \\t]+([0-9]+)$")
    file(STRINGS "${_ZeroMQ_H}" _ZeroMQ_ver REGEX "${_ZeroMQ_expr}")
    string(REGEX MATCH "${_ZeroMQ_expr}" ZeroMQ_ver "${_ZeroMQ_ver}")
    set(${_ZeroMQ_VER_OUTPUT} "${CMAKE_MATCH_1}" PARENT_SCOPE)
endfunction()

_zmqver_EXTRACT("ZMQ_VERSION_MAJOR" ZeroMQ_VERSION_MAJOR)
_zmqver_EXTRACT("ZMQ_VERSION_MINOR" ZeroMQ_VERSION_MINOR)
_zmqver_EXTRACT("ZMQ_VERSION_PATCH" ZeroMQ_VERSION_PATCH)

message(STATUS "ZeroMQ version: ${ZeroMQ_VERSION_MAJOR}.${ZeroMQ_VERSION_MINOR}.${ZeroMQ_VERSION_PATCH}")

# Get library
find_library( ZEROMQ_LIBRARY NAMES zmq libzmq "libzmq-mt-${ZeroMQ_VERSION_MAJOR}_${ZeroMQ_VERSION_MINOR}_${ZeroMQ_VERSION_PATCH}")

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args( ZEROMQ DEFAULT_MSG ZEROMQ_LIBRARY ZEROMQ_INCLUDE_DIR )

mark_as_advanced( ZEROMQ_LIBRARY ZEROMQ_INCLUDE_DIR )
