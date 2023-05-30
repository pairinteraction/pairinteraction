# Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
#
# This file is part of the pairinteraction GUI.
#
# The pairinteraction GUI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The pairinteraction GUI is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the pairinteraction GUI. If not, see <http://www.gnu.org/licenses/>.
import numpy as np
import pyqtgraph as pg


# === Points item (can be used with pyqtgraph) ===


class PointsItem(pg.Qt.QtWidgets.QGraphicsItem):
    def __init__(self, x=None, y=None, size=1, alpha=80, color=(0, 0, 0)):
        pg.Qt.QtWidgets.QGraphicsItem.__init__(self)
        self.size = size
        self.alpha = alpha
        # self.pen = pg.mkPen((0,0,0,self.alpha),width=self.size,style=QtCore.Qt.CustomDashLine)
        # self.pen.setDashPattern([1, 20, 5, 4])
        self.pen = pg.mkPen(color + (self.alpha,), width=self.size, cosmetic=True)
        self.setData(x, y)
        # self.ItemIgnoresTransformations = True
        # self.setFlag(QtWidgets.QGraphicsItem.ItemIgnoresTransformations, True)

    def setData(self, x, y):
        if x is None:
            x = np.array([])
            y = np.array([])
        npoints = len(x)

        # see http://stackoverflow.com/questions/20119777
        self.qpoints = pg.QtGui.QPolygonF(npoints)
        vptr = self.qpoints.data()
        vptr.setsize(np.dtype(float).itemsize * 2 * npoints)
        data = np.ndarray(shape=(npoints, 2), dtype=float, buffer=memoryview(vptr))
        data.setflags(write=True)
        data[:, 0] = x
        data[:, 1] = y

        self.bounds = pg.QtCore.QRectF(x.min(), y.min(), x.max() - x.min(), y.max() - y.min())
        self.prepareGeometryChange()

    def boundingRect(self):
        return self.bounds

    def paint(self, p, *args):
        p.setPen(self.pen)
        p.drawPoints(self.qpoints)


# === Multi line item (can be used with pyqtgraph) ===
# see https://stackoverflow.com/questions/17103698/plotting-large-arrays-in-pyqtgraph


class MultiLine(pg.Qt.QtWidgets.QGraphicsPathItem):
    def __init__(self, x, y, size=1, alpha=80, color=(0, 0, 0)):
        """x and y are 2D arrays of shape (Nplots, Nsamples)"""
        connections = np.ones(x.shape, dtype=bool)
        connections[:, -1] = 0  # don't draw the segment between each trace

        self.path = pg.arrayToQPath(x.flatten(), y.flatten(), connections.flatten())
        pg.Qt.QtWidgets.QGraphicsPathItem.__init__(self, self.path)
        pen = pg.mkPen(color + (alpha,), width=size, cosmetic=True)
        self.setPen(pen)

    # Override because QGraphicsPathItem.shape is too expensive.
    def shape(self):
        return pg.Qt.QtWidgets.QGraphicsItem.shape(self)

    def boundingRect(self):
        return self.path.boundingRect()
