# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
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
    from collections.abc import Sequence

    from numpy.typing import NDArray
    from typing_extensions import Concatenate

    from pairinteraction_gui.calculate.calculate_base import Parameters, Results
    from pairinteraction_gui.page import SimulationPage

logger = logging.getLogger(__name__)


class PlotWidget(WidgetV):
    """Widget for displaying plots with controls."""

    margin = (0, 0, 0, 0)
    spacing = 15

    def __init__(self, parent: SimulationPage) -> None:
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

    parameters: Parameters[Any] | None = None
    results: Results | None = None

    fit_idx: int = 0
    fit_type: str = ""
    fit_data_highlight: mpl.collections.PathCollection
    fit_curve: Sequence[mpl.lines.Line2D]

    def setupWidget(self) -> None:
        super().setupWidget()

        mappable = plt.cm.ScalarMappable(cmap=alphamagma, norm=Normalize(vmin=0, vmax=1))
        self.canvas.fig.colorbar(mappable, ax=self.canvas.ax, label="Overlap with state of interest")
        self.canvas.fig.tight_layout()

    def plot(self, parameters: Parameters[Any], results: Results) -> None:
        self.reset_fit()
        ax = self.canvas.ax
        ax.clear()

        # store data to allow fitting later on
        self.parameters = parameters
        self.results = results

        x_values = parameters.get_x_values()
        energies = results.energies

        try:
            ax.plot(x_values, np.array(energies), c="0.75", lw=0.25, zorder=-10)
        except ValueError as err:
            if "inhomogeneous shape" in str(err):
                for x_value, es in zip(x_values, energies):
                    ax.plot([x_value] * len(es), es, c="0.75", ls="None", marker=".", zorder=-10)
            else:
                raise err

        # Flatten the arrays for scatter plot and repeat x value for each energy
        # (dont use numpy.flatten, etc. to also handle inhomogeneous shapes)
        x_repeated = np.hstack([val * np.ones_like(es) for val, es in zip(x_values, energies)])
        energies_flattened = np.hstack(energies)
        overlaps_flattened = np.hstack(results.ket_overlaps)

        min_overlap = 1e-4
        inds: NDArray[Any] = np.argwhere(overlaps_flattened > min_overlap).flatten()
        inds = inds[np.argsort(overlaps_flattened[inds])]

        if len(inds) > 0:
            ax.scatter(
                x_repeated[inds],
                energies_flattened[inds],
                c=overlaps_flattened[inds],
                s=15,
                vmin=0,
                vmax=1,
                cmap=alphamagma,
            )

        ax.set_xlabel(parameters.get_x_label())
        ax.set_ylabel("Energy [GHz]")

        self.canvas.fig.tight_layout()

    def add_cursor(self, parameters: Parameters[Any], results: Results) -> None:
        # Remove any existing cursors to avoid duplicates
        if hasattr(self, "mpl_cursor"):
            if hasattr(self.mpl_cursor, "remove"):  # type: ignore
                self.mpl_cursor.remove()  # type: ignore
            del self.mpl_cursor  # type: ignore

        energies = results.energies
        x_values = parameters.get_x_values()

        artists = []
        for idx, labels in results.state_labels.items():
            x = x_values[idx]
            for energy, label in zip(energies[idx], labels):
                artist = self.canvas.ax.plot(x, energy, "d", c="0.93", alpha=0.5, ms=7, label=label, zorder=-20)
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

    def fit(self, fit_type: str = "c6") -> None:  # noqa: PLR0912, PLR0915, C901
        """Fits a potential curve and displays the fit values.

        Args:
            fit_type: Type of fit to perform. Options are:
              c6: E = E0 + C6 * r^6
              c3: E = E0 + C3 * r^3
              c3+c6: E = E0 + C3 * r^3 + C6 * r^6

        Iterative calls will iterate through the potential curves

        """
        if self.parameters is None or self.results is None:
            logger.warning("No data to fit.")
            return

        energies = self.results.energies
        x_values = np.array(self.parameters.get_x_values())
        overlaps_list = self.results.ket_overlaps

        fit_func: Callable[Concatenate[NDArray[Any], ...], NDArray[Any]]
        if fit_type == "c6":
            fit_func = fit_c6
            fitlabel = "E0 = {0:.3f} GHz\nC6 = {1:.3f} GHz*µm^6"
        elif fit_type == "c3":
            fit_func = fit_c3
            fitlabel = "E0 = {0:.3f} GHz\nC3 = {1:.3f} GHz*µm^3"
        elif fit_type == "c3+c6":
            fit_func = fit_c3_c6
            fitlabel = "E0 = {0:.3f} GHz\nC3 = {1:.3f} GHz*µm^3\nC6 = {2:.3f} GHz*µm^6"
        else:
            raise ValueError(f"Unknown fit type: {fit_type}")

        # increase the selected potential curve by one if we use the same fit type
        if self.fit_type == fit_type:
            self.fit_idx = (self.fit_idx + 1) % len(energies)
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
        idxs = [np.argpartition(overlaps_list[0], -self.fit_idx)[-self.fit_idx]]
        last_overlap = overlaps_list[0][idxs[-1]]
        for overlaps in overlaps_list[1:]:
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
        energies_fit = np.array([energy[idx] for energy, idx in zip(energies, idxs)])

        # stop highlighting the previous fit
        if hasattr(self, "fit_data_highlight"):
            self.fit_data_highlight.remove()
        if hasattr(self, "fit_curve"):
            for curve in self.fit_curve:
                curve.remove()

        self.fit_data_highlight = self.canvas.ax.scatter(x_values, energies_fit, c="green", s=5)

        try:
            fit_params = curve_fit(fit_func, x_values, energies_fit)[0]
        except (RuntimeError, TypeError):
            logger.warning("Curve fit failed.")
        else:
            self.fit_curve = self.canvas.ax.plot(
                x_values,
                fit_func(x_values, *fit_params),
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


def fit_c3(x: NDArray[Any], /, e0: float, c3: float) -> NDArray[Any]:
    return e0 + c3 / x**3


def fit_c6(x: NDArray[Any], /, e0: float, c6: float) -> NDArray[Any]:
    return e0 + c6 / x**6


def fit_c3_c6(x: NDArray[Any], /, e0: float, c3: float, c6: float) -> NDArray[Any]:
    return e0 + c3 / x**3 + c6 / x**6
