# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QWheelEvent
from PySide6.QtWidgets import QWidget


class MatplotlibCanvas(FigureCanvasQTAgg):
    """Canvas for matplotlib figures."""

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        """Initialize the canvas with a figure."""
        self.fig, self.ax = plt.subplots()
        super().__init__(self.fig)

        self.setup_zoom()

    def setup_zoom(self) -> None:
        """Set up mouse wheel zoom functionality."""
        # Wheel event accumulation variables
        self.wheel_accumulation: float = 0
        self.last_wheel_pos: list[list[float]] = []
        self.wheel_timer = QTimer(self)
        self.wheel_timer.setSingleShot(True)
        self.wheel_timer.timeout.connect(self.apply_accumulated_zoom)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMouseTracking(True)

    def wheelEvent(self, event: QWheelEvent) -> None:
        """Handle mouse wheel events for zooming."""
        self.wheel_accumulation += event.angleDelta().y() / 120
        self.last_wheel_pos.append([event.position().x(), event.position().y()])
        self.wheel_timer.start(100)  # Apply zoom after 100ms of inactivity

    def apply_accumulated_zoom(self) -> None:
        """Apply the accumulated zoom from wheel events."""
        if not self.wheel_accumulation:
            return

        x_min, x_max = self.ax.get_xlim()
        y_min, y_max = self.ax.get_ylim()

        scale_factor = 1 - 0.1 * self.wheel_accumulation

        # Get the mouse position in data coordinates
        wheel_pos_mean = np.mean(self.last_wheel_pos, axis=0)
        x_data, y_data = self.ax.transData.inverted().transform(wheel_pos_mean)
        y_data = -(y_data - (y_max + y_min) / 2) + (y_max + y_min) / 2  # y_data is mirrored bottom / top

        self.wheel_accumulation = 0
        self.last_wheel_pos = []

        if x_data > x_max or y_data > y_max or (x_data < x_min and y_data < y_min):
            return

        if x_min <= x_data <= x_max:
            x_min_new = x_data - (x_data - x_min) * scale_factor
            x_max_new = x_data + (x_max - x_data) * scale_factor
            self.ax.set_xbound(x_min_new, x_max_new)
            self.ax.set_autoscalex_on(False)

        if y_min <= y_data <= y_max:
            y_min_new = y_data - (y_data - y_min) * scale_factor
            y_max_new = y_data + (y_max - y_data) * scale_factor
            self.ax.set_ybound(y_min_new, y_max_new)
            self.ax.set_autoscaley_on(False)

        # Redraw the canvas
        self.draw()
