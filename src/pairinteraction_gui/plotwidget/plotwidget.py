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
from PySide6.QtWidgets import QHBoxLayout

from pairinteraction.visualization.colormaps import alphamagma
from pairinteraction_gui.plotwidget.canvas import MatplotlibCanvas
from pairinteraction_gui.plotwidget.navigation_toolbar import CustomNavigationToolbar
from pairinteraction_gui.qobjects import WidgetV
from pairinteraction_gui.theme import plot_widget_theme

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

        self.setStyleSheet(plot_widget_theme)

    def setupWidget(self) -> None:
        self.canvas = MatplotlibCanvas(self)
        self.navigation_toolbar = CustomNavigationToolbar(self.canvas, self)

        top_layout = QHBoxLayout()
        top_layout.addStretch(1)
        top_layout.addWidget(self.navigation_toolbar)
        self.layout().addLayout(top_layout)

        self.layout().addWidget(self.canvas, stretch=1)

    def clear(self) -> None:
        self.canvas.ax.clear()
        self.canvas.draw()


class PlotEnergies(PlotWidget):
    """Plotwidget for plotting energy levels."""

    def setupWidget(self) -> None:
        super().setupWidget()

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
            ax.plot(x_list, np.array(energies_list), c="0.75", lw=0.25, zorder=-10)
        except ValueError as err:
            if "inhomogeneous shape" in str(err):
                for x_value, es in zip(x_list, energies_list):
                    ax.plot([x_value] * len(es), es, c="0.75", ls="None", marker=".", zorder=-10)
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

        ax.set_xlabel(xlabel)
        ax.set_ylabel("Energy [GHz]")

        self.canvas.fig.tight_layout()

    def add_cursor(
        self, x_value: list[float], energies: list["NDArray[Any]"], state_labels: dict[int, list[str]]
    ) -> None:
        # Remove any existing cursors to avoid duplicates
        if hasattr(self, "mpl_cursor"):
            if hasattr(self.mpl_cursor, "remove"):  # type: ignore
                self.mpl_cursor.remove()  # type: ignore
            del self.mpl_cursor  # type: ignore

        ax = self.canvas.ax

        artists = []
        for idx, labels in state_labels.items():
            x = x_value[idx]
            for energy, label in zip(energies[idx], labels):
                artist = ax.plot(x, energy, "d", c="0.93", alpha=0.5, ms=7, label=label, zorder=-20)
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
            sel.annotation.set_text(label.replace(" + ", "\n + "))
