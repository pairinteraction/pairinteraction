# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
import time
from pathlib import Path
from typing import Any, Optional

import nbformat
from nbconvert import PythonExporter
from PySide6.QtCore import QSize, Qt
from PySide6.QtGui import QHideEvent, QMovie, QShowEvent
from PySide6.QtWidgets import (
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QMenu,
    QPushButton,
    QStyle,
    QToolBox,
)

from pairinteraction_gui.app import Application
from pairinteraction_gui.calculate.calculate_base import Parameters, Results
from pairinteraction_gui.config import BaseConfig
from pairinteraction_gui.plotwidget.plotwidget import PlotEnergies, PlotWidget
from pairinteraction_gui.qobjects import WidgetV
from pairinteraction_gui.qobjects.events import show_status_tip
from pairinteraction_gui.qobjects.item import Item, RangeItem
from pairinteraction_gui.qobjects.named_stacked_widget import NamedStackedWidget
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

    plotwidget: PlotWidget

    _button_style = """
        QPushButton {
            padding: 8px 16px;
            background-color: #343a40;
            color: #ffffff;
            border: none;
            border-radius: 4px;
            font-weight: bold;
            font-size: 14px;
        }
        QPushButton:hover {
            background-color: #495057;
        }
        QPushButton:pressed {
            background-color: #000000;
        }
    """

    _button_menu_style = """
        QMenu {
            background-color: #ffffff;
            border: 1px solid #ffffff;
            border-radius: 4px;
            padding: 4px;
        }
        QMenu::item {
            padding: 6px 24px;
            color: #000000;
            font-size: 14px;
        }
        QMenu::item:selected {
            background-color: #ffffff;
        }
    """

    def setupWidget(self) -> None:
        self.toolbox = QToolBox()

    def postSetupWidget(self) -> None:
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


class CalculationPage(SimulationPage):
    """Base class for all pages with a calculation button."""

    plotwidget: PlotEnergies

    def setupWidget(self) -> None:
        super().setupWidget()

        self.plotwidget = PlotEnergies(self)
        self.layout().addWidget(self.plotwidget)

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

        # Calculate/Abort stacked buttons
        self.calculate_and_abort = NamedStackedWidget[QPushButton](self)

        calculate_button = QPushButton("Calculate")
        calculate_button.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_BrowserReload))
        calculate_button.clicked.connect(self.calculate_clicked)
        calculate_button.setStyleSheet(self._button_style)
        self.calculate_and_abort.addNamedWidget(calculate_button, "Calculate")

        abort_button = QPushButton("Abort")
        abort_button.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_BrowserStop))
        abort_button.clicked.connect(self.abort_clicked)
        abort_button.setStyleSheet(self._button_style)
        self.calculate_and_abort.addNamedWidget(abort_button, "Abort")

        self.calculate_and_abort.setFixedHeight(50)
        control_layout.addWidget(self.calculate_and_abort, stretch=2)

        # Create export button with menu
        export_button = QPushButton("Export")
        export_button.setObjectName("Export")
        export_button.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_DialogSaveButton))
        export_button.setStyleSheet(self._button_style)
        export_menu = QMenu(self)
        export_menu.setStyleSheet(self._button_menu_style)
        export_menu.addAction("Export as PNG", self.export_png)
        export_menu.addAction("Export as Python script", self.export_python)
        export_menu.addAction("Export as Jupyter notebook", self.export_notebook)
        export_button.setMenu(export_menu)
        export_button.setFixedHeight(50)
        control_layout.addWidget(export_button, stretch=1)

        self.layout().addLayout(control_layout)

        self.energy_range = RangeItem(
            self,
            "Calculate the energies from",
            (-999, 999),
            (-0.5, 0.5),
            unit="GHz",
            checked=False,
            tooltip_label="energy",
        )
        self.layout().addWidget(self.energy_range)

        self.fast_mode = Item(self, "Use fast calculation mode", {}, "", checked=True)
        self.layout().addWidget(self.fast_mode)

    def calculate_clicked(self) -> None:
        self.before_calculate()

        def update_plot(
            parameters_and_results: tuple[Parameters[Any], Results],
        ) -> None:
            worker_plot = Worker(self.update_plot, *parameters_and_results)
            worker_plot.start()

        worker = Worker(self.calculate)
        worker.signals.result.connect(update_plot)
        worker.signals.finished.connect(self.after_calculate)
        worker.start()

    def before_calculate(self) -> None:
        show_status_tip(self, "Calculating... Please wait.", logger=logger)
        self.calculate_and_abort.setCurrentNamedWidget("Abort")
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
        else:
            show_status_tip(self, f"Calculation failed after {time_needed:.2f} seconds.", logger=logger)

        self.calculate_and_abort.setCurrentNamedWidget("Calculate")

    def calculate(self) -> tuple[Parameters[Any], Results]:
        raise NotImplementedError("Subclasses must implement this method")

    def update_plot(self, parameters: Parameters[Any], results: Results) -> None:
        energies = results.energies
        overlaps = results.ket_overlaps

        x_values = parameters.get_x_values()
        x_label = parameters.get_x_label()

        self.plotwidget.plot(x_values, energies, overlaps, x_label)

        ind = 0 if parameters.n_atoms == 1 else -1
        self.plotwidget.add_cursor(x_values[ind], energies[ind], results.state_labels_0)

        self.plotwidget.canvas.draw()

    def export_png(self) -> None:
        """Export the current plot as a PNG file."""
        logger.debug("Exporting results as PNG...")

        filename, _ = QFileDialog.getSaveFileName(self, "Save Plot", "", "PNG Files (*.png)")

        if filename:
            filename = filename.removesuffix(".png") + ".png"
            self.plotwidget.canvas.fig.savefig(
                filename, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none"
            )
            logger.info("Plot saved as %s", filename)

    def export_python(self) -> None:
        """Export the current calculation as a Python script."""
        logger.debug("Exporting results as Python script...")

        filename, _ = QFileDialog.getSaveFileName(self, "Save Python Script", "", "Python Files (*.py)")

        if filename:
            filename = filename.removesuffix(".py") + ".py"

            template_path = (
                Path(__file__).parent.parent / "export_templates" / self._get_export_notebook_template_name()
            )
            with Path(template_path).open() as f:
                notebook = nbformat.read(f, as_version=4)

            exporter = PythonExporter(exclude_output_prompt=True, exclude_input_prompt=True)
            content, _ = exporter.from_notebook_node(notebook)

            replacements = self._get_export_replacements()
            for key, value in replacements.items():
                content = content.replace(key, str(value))

            with Path(filename).open("w") as f:
                f.write(content)

            logger.info("Python script saved as %s", filename)

    def export_notebook(self) -> None:
        """Export the current calculation as a Jupyter notebook."""
        logger.debug("Exporting results as Jupyter notebook...")

        filename, _ = QFileDialog.getSaveFileName(self, "Save Jupyter Notebook", "", "Jupyter Notebooks (*.ipynb)")

        if filename:
            filename = filename.removesuffix(".ipynb") + ".ipynb"

            template_path = (
                Path(__file__).parent.parent / "export_templates" / self._get_export_notebook_template_name()
            )
            with Path(template_path).open() as f:
                notebook = nbformat.read(f, as_version=4)

            replacements = self._get_export_replacements()
            for cell in notebook.cells:
                if cell.cell_type == "code":
                    source = cell.source
                    for key, value in replacements.items():
                        source = source.replace(key, str(value))
                    cell.source = source

            nbformat.write(notebook, filename)

            logger.info("Jupyter notebook saved as %s", filename)

    def _get_export_notebook_template_name(self) -> str:
        raise NotImplementedError("Subclasses must implement this method")

    def _get_export_replacements(self) -> dict[str, str]:
        # Override this method in subclasses to provide specific replacements for the export
        return {}

    def abort_clicked(self) -> None:
        """Handle abort button click."""
        logger.debug("Aborting calculation.")
        Application.terminate_all_processes()
        Application.terminate_all_threads()
        self.after_calculate(False)
        show_status_tip(self, "Calculation aborted.", logger=logger)
