import logging
import time
from pathlib import Path
from typing import Optional

from PySide6.QtGui import QHideEvent, QShowEvent
from PySide6.QtWidgets import (
    QHBoxLayout,
    QPushButton,
    QToolBox,
)

from pairinteraction_gui.config import BaseConfig
from pairinteraction_gui.qobjects import WidgetV
from pairinteraction_gui.qobjects.events import show_status_tip
from pairinteraction_gui.worker import THREADPOOL, Worker

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

    def setupWidget(self) -> None:
        self.toolbox = QToolBox()

        # Control panel below the plot
        control_layout = QHBoxLayout()
        calculate_button = QPushButton("Calculate")
        calculate_button.setObjectName("Calculate")
        calculate_button.clicked.connect(self._thread_calculate)
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

    def _thread_calculate(self) -> None:
        self.findChild(QPushButton, "Calculate").setEnabled(False)
        worker = Worker(self.calculate)
        worker.signals.finished.connect(self.calculate_finished)
        THREADPOOL.start(worker)

    def calculate(self) -> None:
        self._start_time = time.perf_counter()
        show_status_tip(self, "Calculating... Please wait.", logger=logger)

    def calculate_finished(self) -> None:
        time_needed = time.perf_counter() - self._start_time
        show_status_tip(self, f"Calculation finished after {time_needed:.2f} seconds.", logger=logger)
        self.findChild(QPushButton, "Calculate").setEnabled(True)
        worker = Worker(self.update_plot)
        THREADPOOL.start(worker)

    def update_plot(self) -> None:
        raise NotImplementedError("Subclasses must implement this method")

    def export(self) -> None:
        logger.debug("Exporting results...")
        return
