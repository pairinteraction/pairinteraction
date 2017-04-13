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

file(READ "${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/plotter.py" TMPTXT)
string(REPLACE "from pyqtgraph import" "from .pyqtgraph import" TMPTXT "${TMPTXT}")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/plotter.py" "${TMPTXT}")

# See https://github.com/pyqtgraph/pyqtgraph/issues/454
file(READ "${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/pyqtgraph/exporters/ImageExporter.py" TMPTXT)
string(REPLACE "np.empty((self.params['width'], self.params['height']," "np.empty((int(self.params['width']), int(self.params['height'])," TMPTXT "${TMPTXT}")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/pairinteraction/pyqtgraph/exporters/ImageExporter.py" "${TMPTXT}")
