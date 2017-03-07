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
# ZMQ_LIBRARY       - the ZeroMQ library
# ZMQ_INCLUDE_DIR   - path including ZeroMQ.h
# ZMQ_BINARY        - ZeroMQ executable

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include (FindPackageHandleStandardArgs)

find_path( ZMQ_INCLUDE_DIR NAMES zmq.h )

find_library( ZMQ_LIBRARY NAMES zmq )

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args( ZMQ DEFAULT_MSG ZMQ_LIBRARY ZMQ_INCLUDE_DIR )

mark_as_advanced( ZMQ_LIBRARY ZMQ_INCLUDE_DIR )
