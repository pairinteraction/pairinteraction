# Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
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

# Find SQLite3
#
# SQLITE3_LIBRARY       - the sqlite3 library
# SQLITE3_INCLUDE_DIR   - path including sqlite3.h
# SQLITE3_BINARY        - sqlite3 executable

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include (FindPackageHandleStandardArgs)

find_path( SQLITE3_INCLUDE_DIR NAMES sqlite3.h )

find_library( SQLITE3_LIBRARY NAMES sqlite3 )

find_program( SQLITE3_BINARY NAMES sqlite3 )

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args( SQLITE3 DEFAULT_MSG SQLITE3_LIBRARY SQLITE3_INCLUDE_DIR SQLITE3_BINARY )

mark_as_advanced( SQLITE3_LIBRARY SQLITE3_INCLUDE_DIR SQLITE3_BINARY )
