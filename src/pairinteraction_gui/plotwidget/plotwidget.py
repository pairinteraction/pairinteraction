# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import mplcursors
import numpy as np
from matplotlib.colors import Normalize
from PySide6.QtWidgets import QPushButton

from pairinteraction.visualization.colormaps import alphamagma
from pairinteraction_gui.plotwidget.canvas import MatplotlibCanvas
from pairinteraction_gui.qobjects import WidgetH, WidgetV
from pairinteraction_gui.qobjects.item import Item, RangeItem

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from pairinteraction_gui.page import SimulationPage

logger = logging.getLogger(__name__)


class PlotWidget(WidgetV):
    """Widget for displaying plots with controls."""

    margin = (0, 0, 0, 0)
    spacing = 15

    def __init__(self, parent: "SimulationPage") -> None:
        """Initialize the base section."""
        mpl.use("Qt5Agg")

        self.page = parent
        super().__init__(parent)

    def setupWidget(self) -> None:
        self.plot_toolbar = WidgetV(self)
        self.canvas = MatplotlibCanvas(self)

    def postSetupWidget(self) -> None:
        top_row = WidgetH(self)

        # Add plot toolbar on left
        top_row.layout().addWidget(self.plot_toolbar)

        # Add reset zoom button on right
        reset_zoom_button = QPushButton("Reset Zoom", self)
        reset_zoom_button.setToolTip(
            "Reset the plot view to its original state. You can zoom in/out using the mousewheel."
        )
        reset_zoom_button.clicked.connect(self.canvas.reset_view)
        top_row.layout().addWidget(reset_zoom_button)

        self.layout().addWidget(top_row)
        self.layout().addWidget(self.canvas, stretch=1)

    def clear(self) -> None:
        self.canvas.ax.clear()
        self.canvas.draw()


class PlotEnergies(PlotWidget):
    """Plotwidget for plotting energy levels."""

    def setupWidget(self) -> None:
        super().setupWidget()

        self.energy_range = RangeItem(
            self,
            "Calculate the energies from",
            (-999, 999),
            (-0.5, 0.5),
            unit="GHz",
            checked=False,
            tooltip_label="energy",
        )
        self.plot_toolbar.layout().addWidget(self.energy_range)

        self.fast_mode = Item(self, "Use fast calculation mode", {}, "", checked=True)
        self.plot_toolbar.layout().addWidget(self.fast_mode)

        mappable = plt.cm.ScalarMappable(cmap=alphamagma, norm=Normalize(vmin=0, vmax=1))
        self.canvas.fig.colorbar(mappable, ax=self.canvas.ax, label="Overlap with state of interest")
        self.canvas.fig.tight_layout()

    def plot(
        self,
        x_list: Sequence[float],
        energies_list: Sequence["NDArray[Any]"],
        overlaps_list: Sequence["NDArray[Any]"],
        xlabel: str,
    ) -> None:
        ax = self.canvas.ax
        ax.clear()

        try:
            ax.plot(x_list, np.array(energies_list), c="0.9", lw=0.25, zorder=-10)
        except ValueError as err:
            if "inhomogeneous shape" in str(err):
                for x_value, es in zip(x_list, energies_list):
                    ax.plot([x_value] * len(es), es, c="0.9", ls="None", marker=".", zorder=-10)
            else:
                raise err

        # Flatten the arrays for scatter plot and repeat x value for each energy
        # (dont use numpy.flatten, etc. to also handle inhomogeneous shapes)
        x_repeated = np.hstack([val * np.ones_like(es) for val, es in zip(x_list, energies_list)])
        energies_flattend = np.hstack(energies_list)
        overlaps_flattend = np.hstack(overlaps_list)

        min_overlap = 1e-4
        inds: NDArray[Any] = np.argwhere(overlaps_flattend > min_overlap).flatten()
        inds = inds[np.argsort(overlaps_flattend[inds])]

        if len(inds) > 0:
            ax.scatter(
                x_repeated[inds],
                energies_flattend[inds],
                c=overlaps_flattend[inds],
                s=15,
                vmin=0,
                vmax=1,
                cmap=alphamagma,
            )

        ylim = ax.get_ylim()
        if abs(ylim[1] - ylim[0]) < 1e-2:
            ax.set_ylim(ylim[0] - 1e-2, ylim[1] + 1e-2)

        ax.set_xlabel(xlabel)
        ax.set_ylabel("Energy [GHz]")
        self.canvas.fig.tight_layout()

    def add_cursor(self, x_value: float, energies: "NDArray[Any]", state_labels_0: list[str]) -> None:
        # Remove any existing cursors to avoid duplicates
        if hasattr(self, "mpl_cursor"):
            if hasattr(self.mpl_cursor, "remove"):  # type: ignore
                self.mpl_cursor.remove()  # type: ignore
            del self.mpl_cursor  # type: ignore

        ax = self.canvas.ax

        artists = []
        for e, ket_label in zip(energies, state_labels_0):
            artist = ax.plot(x_value, e, "o", c="0.9", ms=5, zorder=-20, fillstyle="none", label=ket_label)
            artists.extend(artist)

        self.mpl_cursor = mplcursors.cursor(
            artists,
            hover=False,
            annotation_kwargs={
                "bbox": {"boxstyle": "round,pad=0.5", "fc": "white", "alpha": 0.9, "ec": "gray"},
                "arrowprops": {"arrowstyle": "->", "connectionstyle": "arc3", "color": "gray"},
            },
        )

        @self.mpl_cursor.connect("add")
        def on_add(sel: mplcursors.Selection) -> None:
            label = sel.artist.get_label()
            sel.annotation.set_text(label)
