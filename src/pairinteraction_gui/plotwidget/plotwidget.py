# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Callable

import matplotlib as mpl
import matplotlib.pyplot as plt
import mplcursors
import numpy as np
from matplotlib.colors import Normalize
from PySide6.QtWidgets import QHBoxLayout
from scipy.optimize import curve_fit

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

    x_list: "NDArray[Any]"
    energies_list: Sequence["NDArray[Any]"]
    overlaps_list: Sequence["NDArray[Any]"]
    fit_idx: int = 0
    fit_type: str = ""
    fit_data_highlight: mpl.collections.PathCollection
    fit_curve: Sequence[mpl.lines.Line2D]

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
        # store data to allow fitting later on
        self.x_list = np.array(x_list)
        self.energies_list = energies_list
        self.overlaps_list = overlaps_list
        self.reset_fit()

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

    def fit(self, fit_type: str = "c6") -> None:  # noqa: PLR0912, C901
        """Fits a potential curve and displays the fit values.

        Args:
            fit_type: Type of fit to perform. Options are:
              c6: E = E0 + C6 * r^6
              c3: E = E0 + C3 * r^3
              c3+c6: E = E0 + C3 * r^3 + C6 * r^6

        Iterative calls will iterate through the potential curves

        """
        fit_func: Callable[..., NDArray[Any]]
        if fit_type == "c6":
            fit_func = lambda x, e0, c6: e0 + c6 / x**6  # noqa: E731
            fitlabel = "E0 = {0:.3f} GHz\nC6 = {1:.3f} GHz*µm^6"
        elif fit_type == "c3":
            fit_func = lambda x, e0, c3: e0 + c3 / x**3  # noqa: E731
            fitlabel = "E0 = {0:.3f} GHz\nC3 = {1:.3f} GHz*µm^3"
        elif fit_type == "c3+c6":
            fit_func = lambda x, e0, c3, c6: e0 + c3 / x**3 + c6 / x**6  # noqa: E731
            fitlabel = "E0 = {0:.3f} GHz\nC3 = {1:.3f} GHz*µm^3\nC6 = {2:.3f} GHz*µm^6"
        else:
            raise ValueError(f"Unknown fit type: {fit_type}")

        # first see if we actually have data to fit
        if not (hasattr(self, "x_list") and hasattr(self, "energies_list") and hasattr(self, "overlaps_list")):
            return

        # increase the selected potential curve by one if we use the same fit type
        if self.fit_type == fit_type:
            self.fit_idx = (self.fit_idx + 1) % len(self.energies_list)
        else:
            self.fit_idx = 1
        self.fit_type = fit_type

        # We want to follow the potential curves. The ordering of energies is just by value, so we
        # need to follow the curve somehow. We go right to left, start at the nth largest value, keep our
        # index as long as the difference in overlap is less than a factor 2 or less than 5% total difference.
        # Otherwise, we search until we find an overlap that is less than a factor 2 different.
        # This is of course a simple heuristic, a more sophisticated approach would do some global optimization
        # of the curves. This approach is simple, fast and robust, but curves may e.g. merge.
        # This does not at all take into account the line shapes of the curves. There is also no trade-off
        # between overlap being close and not doing jumps.
        idxs = [np.argpartition(self.overlaps_list[0], -self.fit_idx)[-self.fit_idx]]
        last_overlap = self.overlaps_list[0][idxs[-1]]
        for overlaps in self.overlaps_list[1:]:
            idx = idxs[-1]
            overlap = overlaps[idx]
            if 0.5 * last_overlap < overlap < 2 * last_overlap or abs(overlap - last_overlap) < 0.05:
                # we keep the current index
                idxs.append(idx)
                last_overlap = overlap
            else:
                # we search until we find an overlap that is less than a factor 2 different
                possible_options = np.argwhere(
                    np.logical_and(overlaps > 0.5 * last_overlap, overlaps < 2 * last_overlap)
                ).flatten()
                if len(possible_options) == 0:
                    # there is no state in that range - our best bet is to keep the current index
                    idxs.append(idx)
                    last_overlap = overlap
                else:
                    # we select the closest possible option
                    best_option = np.argmin(np.abs(possible_options - idx))
                    idxs.append(possible_options[best_option])
                    last_overlap = overlaps[idxs[-1]]

        # this could be a call to np.take_along_axis if the sizes match, but the handling of inhomogeneous shapes
        # in the plot() function makes me worry they won't, so I go for a slower python for loop...
        energies = np.array([energy[idx] for energy, idx in zip(self.energies_list, idxs)])

        # stop highlighting the previous fit
        if hasattr(self, "fit_data_highlight"):
            self.fit_data_highlight.remove()
        if hasattr(self, "fit_curve"):
            for curve in self.fit_curve:
                curve.remove()

        self.fit_data_highlight = self.canvas.ax.scatter(self.x_list, energies, c="green", s=5)

        try:
            fit_params = curve_fit(fit_func, self.x_list, energies)[0]
        except (RuntimeError, TypeError):
            logger.warning("Curve fit failed.")
        else:
            self.fit_curve = self.canvas.ax.plot(
                self.x_list,
                fit_func(self.x_list, *fit_params),
                c="green",
                linestyle="dashed",
                lw=2,
                label=fitlabel.format(*fit_params),
            )
            self.canvas.ax.legend()

        self.canvas.draw()

    def clear(self) -> None:
        super().clear()
        self.reset_fit()

    def reset_fit(self) -> None:
        """Clear fit output and reset fit index."""
        # restart at first potential curve
        self.fit_idx = 0
        # and also remove any previous highlighting/fit display
        if hasattr(self, "fit_data_highlight"):
            del self.fit_data_highlight
        if hasattr(self, "fit_curve"):
            del self.fit_curve
