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

file(READ "${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/plotter.py" TMPTXT)
string(REPLACE "from pyqtgraph import" "from .pyqtgraph import" TMPTXT "${TMPTXT}")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/plotter.py" "${TMPTXT}")

# See https://github.com/pyqtgraph/pyqtgraph/issues/454
file(READ "${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/pyqtgraph/exporters/ImageExporter.py" TMPTXT)
string(REPLACE "np.empty((self.params['width'], self.params['height']," "np.empty((int(self.params['width']), int(self.params['height'])," TMPTXT "${TMPTXT}")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/pyqtgraph/exporters/ImageExporter.py" "${TMPTXT}")
