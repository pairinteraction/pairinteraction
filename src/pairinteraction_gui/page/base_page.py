# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
import time
from pathlib import Path
from typing import Any, Optional

from PySide6.QtCore import QSize, Qt
from PySide6.QtGui import QHideEvent, QMovie, QShowEvent
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QPushButton,
    QToolBox,
)

from pairinteraction_gui.calculate.calculate_base import Parameters, Results
from pairinteraction_gui.config import BaseConfig
from pairinteraction_gui.plotwidget.plotwidget import PlotEnergies
from pairinteraction_gui.qobjects import WidgetV
from pairinteraction_gui.qobjects.events import show_status_tip
from pairinteraction_gui.worker import Worker

logger = logging.getLogger(__name__)


class BasePage(WidgetV):
    """Base class for all pages in this application."""

    margin = (20, 20, 20, 20)
    spacing = 15

    title: str
    tooltip: str
    icon_path: Optional[Path] = None

    def showEvent(self, event: QShowEvent) -> None:
        """Show event."""
        super().showEvent(event)
        self.window().setWindowTitle("Pairinteraction - " + self.title)


class SimulationPage(BasePage):
    """Base class for all simulation pages in this application."""

    parameters: Parameters
    results: Results
    plotwidget: PlotEnergies

    def setupWidget(self) -> None:
        self.toolbox = QToolBox()

        # Setup loading animation
        self.loading_label = QLabel(self)
        gif_path = Path(__file__).parent.parent / "images" / "loading.gif"
        self.loading_movie = QMovie(str(gif_path))
        self.loading_movie.setScaledSize(QSize(100, 100))  # Make the gif larger
        self.loading_label.setMovie(self.loading_movie)
        self.loading_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.loading_label.hide()

        # Control panel below the plot
        control_layout = QHBoxLayout()
        calculate_button = QPushButton("Calculate")
        calculate_button.setObjectName("Calculate")
        calculate_button.clicked.connect(self.calculate_clicked)
        control_layout.addWidget(calculate_button)

        # export_button = QPushButton("Export")
        # export_button.setObjectName("Export")
        # export_button.clicked.connect(self.export)
        # control_layout.addWidget(export_button)

        self.layout().addLayout(control_layout)

    def postSetupWidget(self) -> None:
        self.layout().addStretch()
        for attr in self.__dict__.values():
            if isinstance(attr, BaseConfig):
                self.toolbox.addItem(attr, attr.title)

    def showEvent(self, event: QShowEvent) -> None:
        super().showEvent(event)
        self.window().dockwidget.setWidget(self.toolbox)
        self.window().dockwidget.setVisible(True)
        self.toolbox.show()

    def hideEvent(self, event: QHideEvent) -> None:
        super().hideEvent(event)
        self.window().dockwidget.setVisible(False)

    def calculate_clicked(self) -> None:
        self.before_calculate()

        worker = Worker(self.calculate)
        worker.signals.finished.connect(self.after_calculate)
        worker.start()

    def before_calculate(self) -> None:
        show_status_tip(self, "Calculating... Please wait.", logger=logger)
        self.findChild(QPushButton, "Calculate").setEnabled(False)
        self.plotwidget.clear()

        # run loading gif
        self.loading_label.setGeometry((self.width() - 100) // 2, (self.height() - 100) // 2, 100, 100)
        self.loading_label.show()
        self.loading_movie.start()

        self._start_time = time.perf_counter()

    def after_calculate(self, success: bool) -> None:
        time_needed = time.perf_counter() - self._start_time

        # stop loading gif
        self.loading_movie.stop()
        self.loading_label.hide()

        if success:
            show_status_tip(self, f"Calculation finished after {time_needed:.2f} seconds.", logger=logger)
            worker = Worker(self.update_plot)
            worker.start()
        else:
            show_status_tip(self, f"Calculation failed after {time_needed:.2f} seconds.", logger=logger)

        self.findChild(QPushButton, "Calculate").setEnabled(True)

    def calculate(self) -> None:
        raise NotImplementedError("Subclasses must implement this method")

    def update_plot(self) -> None:
        energies = self.results.energies
        overlaps = self.results.ket_overlaps

        x_values = self.parameters.get_x_values()
        x_label = self.parameters.get_x_label()

        self.plotwidget.plot(x_values, energies, overlaps, x_label)

        ind = 0 if self.parameters.n_atoms == 1 else -1
        self.plotwidget.add_cursor(x_values[ind], energies[ind], self.results.state_labels_0)

        self.plotwidget.canvas.draw()

    def export(self) -> None:
        raise NotImplementedError("Subclasses must implement this method")
