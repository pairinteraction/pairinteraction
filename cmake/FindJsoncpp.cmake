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

# Find JsonCPP
#
# JSONCPP_LIBRARIES      - the jsoncpp library
# JSONCPP_INCLUDE_DIRS   - directory with all header files

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include (FindPackageHandleStandardArgs)

find_path( JSONCPP_INCLUDE_DIRS NAMES json/json.h PATH_SUFFIXES jsoncpp )

find_library( JSONCPP_LIBRARIES NAMES jsoncpp )

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args( JSONCPP DEFAULT_MSG JSONCPP_LIBRARIES JSONCPP_INCLUDE_DIRS )

mark_as_advanced( JSONCPP_LIBRARIES JSONCPP_INCLUDE_DIRS )
