import logging
from collections.abc import Sequence
from itertools import chain
from typing import TYPE_CHECKING, Any

import matplotlib
import numpy as np
from PySide6.QtWidgets import QPushButton

from pairinteraction_gui.plotwidget.canvas import MatplotlibCanvas
from pairinteraction_gui.qobjects import WidgetH, WidgetV
from pairinteraction_gui.qobjects.item import Item, RangeItem
from pairinteraction_gui.qobjects.spin_boxes import DoubleSpinBox

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from pairinteraction_gui.config.system import RangeObject
    from pairinteraction_gui.page.base_page import SimulationPage

logger = logging.getLogger(__name__)
matplotlib.use("Qt5Agg")


class PlotWidget(WidgetV):
    """Widget for displaying plots with controls."""

    margin = (0, 0, 0, 0)
    spacing = 15

    def __init__(self, parent: "SimulationPage") -> None:
        """Initialize the base section."""
        self.page = parent
        super().__init__(parent)

    def setupWidget(self) -> None:
        self.plot_toolbar = WidgetV(self)
        self.canvas = MatplotlibCanvas(self)

    def postSetupWidget(self) -> None:
        top_row = WidgetH(self)

        # Add plot toolbar on left
        if self.plot_toolbar.layout().count():
            self.plot_toolbar.layout().addStretch()
        top_row.layout().addWidget(self.plot_toolbar)

        # Add reset zoom button on right
        reset_zoom_button = QPushButton("Reset Zoom", self)
        reset_zoom_button.setToolTip(
            "Reset the plot view to its original state. You can zoom in/out using the mousewheel."
        )
        reset_zoom_button.clicked.connect(self.canvas.reset_view)
        top_row.layout().addWidget(reset_zoom_button)

        self.layout().addWidget(top_row)
        self.layout().addWidget(self.canvas)

    def clear(self) -> None:
        self.canvas.ax.clear()
        self.canvas.draw()


class PlotEnergies(PlotWidget):
    """Plotwidget for plotting energy levels."""

    def setupWidget(self) -> None:
        super().setupWidget()

        min_spinbox = DoubleSpinBox(self, vmin=-999, vmax=999, vdefault=-0.5, decimals=2)
        max_spinbox = DoubleSpinBox(self, vmin=-999, vmax=999, vdefault=0.5, decimals=2)
        self.energy_range = RangeItem(
            self, "Calculate the energies from", min_spinbox, max_spinbox, "GHz", checked=False
        )
        self.plot_toolbar.layout().addWidget(self.energy_range)

        self.fast_mode = Item(self, "Use fast calculation mode", {}, "", checked=False)
        self.plot_toolbar.layout().addWidget(self.fast_mode)

    def plot(
        self,
        energies: Sequence["NDArray[Any]"],
        x_ranges: dict[str, "RangeObject"],
        overlaps: Sequence["NDArray[Any]"],
    ) -> None:
        ax = self.canvas.ax
        ax.clear()

        x_values, xlabel = self._get_x_values_and_label_from_ranges(x_ranges)

        try:
            ax.plot(x_values, np.array(energies), c="0.9", lw=0.25, zorder=-10)
        except ValueError as err:
            if "inhomogeneous shape" in str(err):
                for x, es in zip(x_values, energies):
                    ax.plot([x] * len(es), es, c="0.9", ls="None", marker=".", zorder=-10)
            else:
                raise err

        _x = chain.from_iterable([x] * len(es) for x, es in zip(x_values, energies))
        x = np.array(list(_x))
        _y = chain.from_iterable(energies)
        y = np.array(list(_y))
        _o = chain.from_iterable(overlaps)
        o = np.array(list(_o))

        min_overlap = 0.0001
        inds: NDArray[Any] = np.argwhere(o > min_overlap).flatten()
        inds = inds[np.argsort(o[inds])]

        if len(inds) > 0:
            log_o = np.log(o[inds])
            alpha: NDArray[Any]
            if log_o.max() - log_o.min() < 1e-10:
                alpha = np.ones_like(log_o)
            else:
                alpha = 1 - log_o / np.log(min_overlap)
                alpha[alpha < 0] = 0
                alpha[alpha > 1] = 1

            ax.scatter(x[inds], y[inds], c=o[inds], alpha=alpha, s=15, vmin=0, vmax=1, cmap="magma_r")

        ylim = ax.get_ylim()
        if abs(ylim[1] - ylim[0]) < 1e-2:
            ax.set_ylim(ylim[0] - 1e-2, ylim[1] + 1e-2)

        ax.set_xlabel(xlabel)
        ax.set_ylabel("Energy [GHz]")

    @staticmethod
    def _get_x_values_and_label_from_ranges(x_ranges: dict[str, "RangeObject"]) -> tuple["NDArray[Any]", str]:
        range_diffs = {key: abs(r[-1] - r[0]) for key, r in x_ranges.items()}
        max_key = max(range_diffs, key=range_diffs.get)  # type: ignore
        x = x_ranges[max_key].list

        xlabel = max_key
        x_units = {"E": "V/cm", "B": "Gauss", "distance": r"$\mu$m", "angle": r"$^\circ$"}
        xlabel += f" [{next((unit for key, unit in x_units.items() if max_key.startswith(key)), '')}]"

        non_constant_keys = [k for k, v in range_diffs.items() if k != max_key and v != 0]
        if non_constant_keys:
            xlabel += f"  ({', '.join(non_constant_keys)} did also change)"

        return x, xlabel
