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

# Find PyUIC
#
# PYUIC_BINARY       - the pyuic GUI builder

# Include these modules to handle the QUIETLY and REQUIRED arguments.
include (FindPackageHandleStandardArgs)

find_program( PYUIC_BINARY NAMES py3uic5 pyuic5 py3uic5.bat pyuic5.bat )

# Set the FOUND variable to TRUE if all listed variables are set.
find_package_handle_standard_args( PyUIC DEFAULT_MSG PYUIC_BINARY )

mark_as_advanced( PYUIC_BINARY )
