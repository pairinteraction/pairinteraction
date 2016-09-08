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

from . import pyqtgraph as pg
import numpy as np


# === Points item (can be used with pyqtgraph) ===

class PointsItem(pg.QtGui.QGraphicsItem):

    def __init__(self, x=None, y=None, size=1, alpha=80, color=(0, 0, 0)):
        pg.QtGui.QGraphicsItem.__init__(self)
        self.size = size
        self.alpha = alpha
        # self.pen = pg.mkPen((0,0,0,self.alpha),width=self.size,style=QtCore.Qt.CustomDashLine)
        # self.pen.setDashPattern([1, 20, 5, 4])
        self.pen = pg.mkPen(color + (self.alpha,),
                            width=self.size, cosmetic=True)
        self.setData(x, y)
        # self.ItemIgnoresTransformations = True
        # self.setFlag(QtGui.QGraphicsItem.ItemIgnoresTransformations, True)

    def setData(self, x, y):
        if x is None:
            x = np.array([])
            y = np.array([])
        npoints = len(x)

        # see http://stackoverflow.com/questions/20119777
        self.qpoints = pg.QtGui.QPolygonF(npoints)
        vptr = self.qpoints.data()
        vptr.setsize(np.dtype(np.float).itemsize * 2 * npoints)
        data = np.ndarray(shape=(npoints, 2), dtype=np.float,
                          buffer=memoryview(vptr))
        data.setflags(write=True)
        data[:, 0] = x
        data[:, 1] = y

        self.bounds = pg.QtCore.QRectF(
            x.min(), y.min(), x.max() - x.min(), y.max() - y.min())
        self.prepareGeometryChange()

    def boundingRect(self):
        return self.bounds

    def paint(self, p, *args):
        p.setPen(self.pen)
        p.drawPoints(self.qpoints)


# === Multi line item (can be used with pyqtgraph) ===
# see https://stackoverflow.com/questions/17103698/plotting-large-arrays-in-pyqtgraph

class MultiLine(pg.QtGui.QGraphicsPathItem):

    def __init__(self, x, y, size=1, alpha=80, color=(0, 0, 0)):
        """x and y are 2D arrays of shape (Nplots, Nsamples)"""
        connections = np.ones(x.shape, dtype=bool)
        connections[:, -1] = 0  # don't draw the segment between each trace

        self.path = pg.arrayToQPath(
            x.flatten(), y.flatten(), connections.flatten())
        pg.QtGui.QGraphicsPathItem.__init__(self, self.path)
        pen = pg.mkPen(color + (alpha,), width=size, cosmetic=True)
        self.setPen(pen)

    # Override because QGraphicsPathItem.shape is too expensive.
    def shape(self):
        return pg.QtGui.QGraphicsItem.shape(self)

    def boundingRect(self):
        return self.path.boundingRect()