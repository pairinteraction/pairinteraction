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

find_program(SQLITE3_BINARY NAMES sqlite3 sqlite REQUIRED)

add_custom_target(database ALL)
add_custom_command(TARGET database
  COMMAND ${SQLITE3_BINARY} ${CMAKE_CURRENT_BINARY_DIR}/quantum_defects.db < ${CMAKE_CURRENT_SOURCE_DIR}/quantum_defects.sql
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/testsuite/pairinteraction/databases
  COMMAND ${SQLITE3_BINARY} ${CMAKE_BINARY_DIR}/testsuite/pairinteraction/databases/quantum_defects.db < ${CMAKE_CURRENT_SOURCE_DIR}/quantum_defects.sql
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/pairinteraction.jl/test/pairinteraction/databases
  COMMAND ${SQLITE3_BINARY} ${CMAKE_BINARY_DIR}/pairinteraction.jl/test/pairinteraction/databases/quantum_defects.db < ${CMAKE_CURRENT_SOURCE_DIR}/quantum_defects.sql
)
add_dependencies(check database)

if( WITH_DMG )

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/quantum_defects.db DESTINATION pairinteraction/databases)

elseif( SKBUILD )

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/quantum_defects.db DESTINATION ${SKBUILD_PLATLIB_DIR}/pairinteraction)

else( )

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/quantum_defects.db DESTINATION share/pairinteraction/pairinteraction/databases)

endif( )
