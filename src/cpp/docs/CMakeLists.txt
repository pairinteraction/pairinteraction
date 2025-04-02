# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

find_package(Doxygen)

if(DOXYGEN_FOUND)
  set(DOXYGEN_CMD SOURCE_DIR=${CMAKE_SOURCE_DIR} ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
  set(DOXYGEN_COMMENT "Building HTML documentation with Doxygen")
else()
  set(DOXYGEN_CMD ${CMAKE_COMMAND} -E echo "Doxygen not found. Documentation will not be generated.")
  set(DOXYGEN_COMMENT "Doxygen not found. Documentation will not be generated.")
endif()

add_custom_target(
  doxygen
  COMMAND ${DOXYGEN_CMD}
  COMMENT "${DOXYGEN_COMMENT}")
