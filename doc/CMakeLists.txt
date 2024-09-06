# SPDX-FileCopyrightText: 2024 Sebastian Weber, Henri Menke, and contributors. All rights reserved.
# SPDX-License-Identifier: GPL-3.0-or-later OR LGPL-3.0-or-later

add_subdirectory(doxygen)
# add_subdirectory(sphinx)  # TODO remove this in new PR

add_custom_target(doc DEPENDS doxygen sphinx)
