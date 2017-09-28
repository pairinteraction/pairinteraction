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

find_package(Sqlite3 REQUIRED)

add_custom_target(database ALL)
add_custom_command(TARGET database
  COMMAND ${SQLITE3_BINARY} ${CMAKE_CURRENT_BINARY_DIR}/quantum_defects.db < ${CMAKE_CURRENT_SOURCE_DIR}/quantum_defects.sql
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/testsuite/pairinteraction/databases
  COMMAND ${SQLITE3_BINARY} ${CMAKE_BINARY_DIR}/testsuite/pairinteraction/databases/quantum_defects.db < ${CMAKE_CURRENT_SOURCE_DIR}/quantum_defects.sql
)

if( WITH_DMG )

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/quantum_defects.db DESTINATION pairinteraction/databases)

else( )

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/quantum_defects.db DESTINATION share/pairinteraction/pairinteraction/databases)

endif( )
