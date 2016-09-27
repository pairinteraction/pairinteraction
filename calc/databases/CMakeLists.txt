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

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(CMAKE_MACOSX_RPATH 0)

find_package(Sqlite3 REQUIRED)
if(SQLITE3_FOUND)
else(SQLITE3_FOUND)
  message(FATAL_ERROR "==================\n"
                      "Sqlite3 not found!\n"
                      "==================\n")
endif(SQLITE3_FOUND)

add_custom_target(database ALL)
add_custom_command(TARGET database
  COMMAND ${SQLITE3_BINARY} ${CMAKE_CURRENT_BINARY_DIR}/quantum_defects.db < ${CMAKE_CURRENT_SOURCE_DIR}/quantum_defects.sql
)

if(DMG)

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/quantum_defects.db DESTINATION databases)

else(DMG)

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/quantum_defects.db DESTINATION share/pairinteraction/calc/databases)

endif(DMG)