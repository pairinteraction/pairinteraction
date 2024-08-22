# SPDX-FileCopyrightText: 2024 Sebastian Weber, Henri Menke, and contributors. All rights reserved.
# SPDX-License-Identifier: GPL-3.0-or-later OR LGPL-3.0-or-later

find_package(Doxygen)

add_custom_target(
  doxygen
  COMMAND SOURCE_DIR=${CMAKE_SOURCE_DIR} ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
  COMMENT "Building HTML documentation with Doxygen")
