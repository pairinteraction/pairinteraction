cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(CMAKE_MACOSX_RPATH 0)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../modules)

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

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/quantum_defects.db DESTINATION share/pairinteraction/)
