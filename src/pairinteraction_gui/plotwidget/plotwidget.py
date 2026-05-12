# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import contextlib
import logging
from typing import TYPE_CHECKING, Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from PySide6.QtGui import QPalette
from PySide6.QtWidgets import QHBoxLayout
from scipy.optimize import curve_fit

from pairinteraction.state.state_atom import StateAtom
from pairinteraction.visualization.colormaps import alphamagma
from pairinteraction_gui.plotwidget.canvas import MatplotlibCanvas
from pairinteraction_gui.plotwidget.navigation_toolbar import CustomNavigationToolbar
from pairinteraction_gui.qobjects import WidgetV
from pairinteraction_gui.qobjects.events import show_status_tip
from pairinteraction_gui.theme import theme_manager

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence
    from typing import Concatenate

    from numpy.typing import NDArray

    from pairinteraction.state import StateBase
    from pairinteraction_gui.calculate.calculate_base import Parameters, Results
    from pairinteraction_gui.calculate.calculate_c6 import ParametersC6, ResultsC6
    from pairinteraction_gui.calculate.calculate_lifetimes import KetData, ParametersLifetimes, ResultsLifetimes
    from pairinteraction_gui.page import SimulationPage

logger = logging.getLogger(__name__)


class PlotWidget(WidgetV):
    """Widget for displaying plots with controls."""

    margin = (0, 0, 0, 0)
    spacing = 15
    _annotations: dict[Any, mpl.text.Annotation]
    parameters: Any | None
    results: Any | None

    def __init__(self, parent: SimulationPage) -> None:
        """Initialize the base section."""
        mpl.use("Qt5Agg")
        self.page = parent
        super().__init__(parent)

        self._annotations = {}
        self._click_cid: int | None = None

    def setupWidget(self) -> None:
        self.canvas = MatplotlibCanvas(self)
        self.navigation_toolbar = CustomNavigationToolbar(self.canvas, self)
        self.navigation_toolbar.setObjectName("PlotNavigationToolBar")

        top_layout = QHBoxLayout()
        top_layout.addStretch(1)
        top_layout.addWidget(self.navigation_toolbar)
        self.layout().addLayout(top_layout)

        self.layout().addWidget(self.canvas, stretch=1)

    def plot(self, parameters: Any, results: Any) -> None:
        self.clear()

        show_status_tip(self, "Plotting...")
        self.canvas.ax.set_xmargin(0)

        self.results = results
        self.parameters = parameters

    def setup_annotations(self, parameters: Any, results: Any) -> None:
        """Add click-based annotations to the plot."""
        show_status_tip(self, "Adding annotations...")

        self.disconnect_click()
        self._click_cid = self.canvas.mpl_connect("button_press_event", self.on_click)  # type: ignore [arg-type]

        if self.clear_annotations not in self.navigation_toolbar._home_callbacks:
            self.navigation_toolbar._home_callbacks += [self.clear_annotations]

    def on_click(self, event: mpl.backend_bases.MouseEvent) -> None:
        """Handle click events on the plot."""
        raise NotImplementedError("Subclasses must implement on_click()")

    def clear(self) -> None:
        self.canvas.ax.clear()
        self.canvas.draw_idle()
        self.clear_annotations()
        self.disconnect_click()
        self.results = None
        self.parameters = None

    def clear_annotations(self) -> None:
        for ann in self._annotations.values():
            with contextlib.suppress(NotImplementedError):
                ann.remove()  # artist may already be gone if ax.clear() was called
        self._annotations.clear()
        self.canvas.draw_idle()

    def disconnect_click(self) -> None:
        if self._click_cid is not None:
            self.canvas.mpl_disconnect(self._click_cid)
            self._click_cid = None


class PlotEnergies(PlotWidget):
    """Plotwidget for plotting energy levels."""

    _annotations: dict[int, mpl.text.Annotation]
    parameters: Parameters[Any] | None
    results: Results | None

    def __init__(self, parent: SimulationPage) -> None:
        super().__init__(parent)

        self.fit_idx = 0
        self.fit_type = ""
        self.fit_data_highlight: mpl.collections.PathCollection | None = None
        self.fit_curve: Sequence[mpl.lines.Line2D] | None = None

    def setupWidget(self) -> None:
        super().setupWidget()

        window_color = theme_manager.get_palette().color(QPalette.ColorRole.Window).name()

        self.canvas.fig.set_facecolor(window_color)
        self.canvas.fig.set_layout_engine(
            "constrained",
            w_pad=0.2,
            h_pad=0.2,
            wspace=0.0,
            hspace=0.0,
        )
        mappable = plt.cm.ScalarMappable(cmap=alphamagma, norm=Normalize(vmin=0, vmax=1))
        cbar = self.canvas.fig.colorbar(mappable, ax=self.canvas.ax, label="Overlap with state of interest", aspect=60)
        cbar.ax.set_zorder(0)
        self.canvas.ax.set_zorder(1)

    def plot(self, parameters: Parameters[Any], results: Results) -> None:
        super().plot(parameters, results)

        x_values = parameters.get_x_values()
        energies = results.energies

        ax = self.canvas.ax
        if len({len(es) for es in energies}) <= 1:  # check if homogeneous shape
            ax.plot(x_values, np.array(energies), c="0.75", lw=0.25, zorder=-10)
        else:
            for x_value, es in zip(x_values, energies, strict=True):
                ax.plot([x_value] * len(es), es, c="0.75", ls="None", marker=".", zorder=-10)

        show_status_tip(self, "Plotting overlaps...")

        # Flatten the arrays for scatter plot and repeat x value for each energy
        # (dont use numpy.flatten, etc. to also handle inhomogeneous shapes)
        x_repeated = np.hstack([val * np.ones_like(es) for val, es in zip(x_values, energies, strict=True)])
        energies_flattened = np.hstack(energies)
        overlaps_flattened = np.hstack(results.ket_overlaps)

        min_overlap = 1e-4
        inds: NDArray[Any] = np.argwhere(overlaps_flattened > min_overlap).ravel()
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
        ax.set_ylabel("Energy (GHz)")

    def setup_annotations(self, parameters: Parameters[Any], results: Results) -> None:
        """Connect click-based state annotation to the energy plot."""
        super().setup_annotations(parameters, results)

        energies = results.energies
        overlaps = results.ket_overlaps
        x_values = parameters.get_x_values()

        self._point_index_map: list[tuple[int, int]] = []
        all_x: list[float] = []
        all_y: list[float] = []
        all_overlaps: list[float] = []
        for idx in range(len(energies)):
            x = x_values[idx]
            for idstate, (energy, overlap) in enumerate(zip(energies[idx], overlaps[idx], strict=True)):
                all_x.append(x)
                all_y.append(float(energy))
                all_overlaps.append(float(overlap))
                self._point_index_map.append((idx, idstate))
        self.pts_data = np.column_stack([all_x, all_y]) if all_x else np.empty((0, 2))
        self.pts_overlaps = np.array(all_overlaps)

    def on_click(self, event: mpl.backend_bases.MouseEvent) -> None:
        if (
            event.inaxes is not self.canvas.ax
            or event.button not in [1, 3]
            or len(self.pts_data) == 0
            or self.navigation_toolbar.mode
        ):
            return
        if event.button == 3:  # right click clears annotations
            self.clear_annotations()
            return

        pts_pos = self.canvas.ax.transData.transform(self.pts_data)
        click_pos = np.array([event.x, event.y])
        dists = np.hypot(pts_pos[:, 0] - click_pos[0], pts_pos[:, 1] - click_pos[1])
        candidates = np.flatnonzero(dists <= 10)  # threshold in pixels
        if len(candidates) == 0:
            self.clear_annotations()
            return
        selected = int(candidates[np.argmax(self.pts_overlaps[candidates])])
        if selected in self._annotations:
            self._annotations[selected].remove()
            del self._annotations[selected]
            self.canvas.draw_idle()
            return
        if self.results is None:
            return
        idstep, idstate = self._point_index_map[selected]
        state: StateBase[Any] = self.results.systems[idstep].get_eigenbasis().get_state(idstate)
        label = state.get_label().replace(" + ", "\n + ").replace(" - ", "\n - ")
        xlim = self.canvas.ax.get_xlim()
        ylim = self.canvas.ax.get_ylim()
        x_frac = (self.pts_data[selected, 0] - xlim[0]) / (xlim[1] - xlim[0])
        y_frac = (self.pts_data[selected, 1] - ylim[0]) / (ylim[1] - ylim[0])
        x_offset = -100 if isinstance(state, StateAtom) else -250
        x_offset = x_offset if x_frac > 0.5 else 0
        y_offset = 15 + 10 * label.count("\n")
        y_offset = -y_offset if y_frac > 0.5 else y_offset
        ann = self.canvas.ax.annotate(
            label,
            xy=(self.pts_data[selected, 0], self.pts_data[selected, 1]),
            xytext=(x_offset, y_offset),
            textcoords="offset points",
            va="center",
            bbox={"boxstyle": "round,pad=0.5", "fc": "white", "alpha": 0.9, "ec": "gray"},
            arrowprops={"arrowstyle": "->", "connectionstyle": "arc3", "color": "gray"},
            clip_on=False,
        )
        ann.set_in_layout(False)
        self._annotations[selected] = ann
        self.canvas.draw_idle()

    def get_annotation_text(self, key: Any) -> str:
        """Return the annotation text for the given bar key."""
        raise NotImplementedError("Subclasses must implement get_annotation_text()")

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
        energies_fit = np.array([energy[idx] for energy, idx in zip(energies, idxs, strict=True)])

        # stop highlighting the previous fit
        if self.fit_data_highlight is not None:
            self.fit_data_highlight.remove()
        if self.fit_curve is not None:
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

        self.canvas.draw_idle()

    def clear(self) -> None:
        super().clear()
        self.reset_fit()

    def reset_fit(self) -> None:
        """Clear fit output and reset fit index."""
        # restart at first potential curve
        self.fit_idx = 0
        # and also remove any previous highlighting/fit display
        self.fit_data_highlight = None
        self.fit_curve = None


class PlotBars(PlotWidget):
    """Plotwidget for plotting bar charts."""

    _annotations: dict[Any, mpl.text.Annotation]
    _bar_data: dict[Any, mpl.patches.Rectangle]

    def setupWidget(self) -> None:
        super().setupWidget()

        window_color = theme_manager.get_palette().color(QPalette.ColorRole.Window).name()
        self.canvas.fig.set_facecolor(window_color)
        self.canvas.fig.set_layout_engine(
            "constrained",
            w_pad=0.2,
            h_pad=0.2,
            wspace=0.0,
            hspace=0.0,
        )

    def on_click(self, event: mpl.backend_bases.MouseEvent) -> None:
        if event.inaxes is not self.canvas.ax or event.button not in [1, 3] or self.navigation_toolbar.mode:
            return
        if event.button == 3:  # right click clears annotations
            self.clear_annotations()
            return

        if event.xdata is None or event.ydata is None:
            return
        click_coords = np.array([event.xdata, event.ydata])

        for _key, rect in self._bar_data.items():
            x0, x1 = rect.get_x(), rect.get_x() + rect.get_width()
            y0, y1 = rect.get_y(), rect.get_y() + rect.get_height()
            inside_x = min(x0, x1) <= click_coords[0] <= max(x0, x1)
            inside_y = min(y0, y1) <= click_coords[1] <= max(y0, y1)
            if inside_x and inside_y:
                key = _key
                break
        else:  # no break -> no bar found
            self.clear_annotations()
            return

        # if we click the same bar again, remove the annotation
        if key in self._annotations:
            self._annotations[key].remove()
            del self._annotations[key]
            self.canvas.draw_idle()
            return

        text = self.get_annotation_text(key)

        xlim = self.canvas.ax.get_xlim()
        ylim = self.canvas.ax.get_ylim()
        bar_cx = rect.get_x() + rect.get_width() / 2
        bar_top = rect.get_y() + rect.get_height()
        x_frac = (bar_cx - xlim[0]) / (xlim[1] - xlim[0])
        y_frac = (bar_top - ylim[0]) / (ylim[1] - ylim[0])
        x_offset = -200 if x_frac > 0.5 else 25
        y_offset = -50 if y_frac > 0.5 else 50
        ann = self.canvas.ax.annotate(
            text,
            xy=(bar_cx, bar_top),
            xytext=(x_offset, y_offset),
            textcoords="offset points",
            bbox={"boxstyle": "round,pad=0.5", "fc": "white", "alpha": 0.9, "ec": "gray"},
            arrowprops={"arrowstyle": "->", "connectionstyle": "arc3", "color": "gray"},
            clip_on=False,
        )

        ann.set_in_layout(False)
        self._annotations[key] = ann
        self.canvas.draw_idle()

    def get_annotation_text(self, key: Any) -> str:
        """Return the annotation text for the given bar key."""
        raise NotImplementedError("Subclasses must implement get_annotation_text()")


class PlotLifetimes(PlotBars):
    """Plotwidget for plotting lifetime/transition rate bar charts."""

    _annotations: dict[tuple[str, int], mpl.text.Annotation]
    _bar_data: dict[tuple[str, int], mpl.patches.Rectangle]
    parameters: ParametersLifetimes | None
    results: ResultsLifetimes | None

    def plot(self, parameters: ParametersLifetimes, results: ResultsLifetimes) -> None:
        super().plot(parameters, results)

        ax = self.canvas.ax

        show_status_tip(self, "Preparing transition rates...")
        labels = ["Spontaneous Decay", "Black Body Radiation"]
        n_list = np.arange(0, np.max([s.n for s in results.kets_bbr + results.kets_sp] + [0]) + 1)
        sorted_rates: dict[str, dict[int, list[tuple[KetData, float]]]] = {}
        for key, kets, rates in [
            (labels[0], results.kets_sp, results.transition_rates_sp),
            (labels[1], results.kets_bbr, results.transition_rates_bbr),
        ]:
            sorted_rates[key] = {n: [] for n in n_list}
            for i, s in enumerate(kets):
                sorted_rates[key][s.n].append((s, rates[i]))
        self.sorted_rates = sorted_rates
        rates_summed = {key: [sum(r for _, r in sorted_rates[key][n]) for n in n_list] for key in sorted_rates}

        show_status_tip(self, "Plotting transition rates...")
        self._bar_data = {}
        for label, color in zip(labels, ["blue", "red"], strict=True):
            artist = ax.bar(n_list, rates_summed[label], label=label, color=color, alpha=0.8)
            for n, rect in zip(n_list, artist.patches, strict=True):
                self._bar_data[(label, n)] = rect

        ax.legend()

        ax.set_xlabel("Principal Quantum Number $n$")
        ax.set_ylabel(r"Transition Rates (1 / ms)")

    def get_annotation_text(self, key: tuple[str, int]) -> str:
        label, n = key
        state_text = "\n".join(f"  - {s}: {r:.5f}/ms" for (s, r) in self.sorted_rates[label][n])
        return f"{label} to n={n}:\n{state_text}"


class PlotC6(PlotBars):
    """Plotwidget for plotting C6 contributions binned by energy gap."""

    BIN_WIDTH_FRACTION = 0.05

    _annotations: dict[int, mpl.text.Annotation]
    _bar_data: dict[int, mpl.patches.Rectangle]
    parameters: ParametersC6 | None
    results: ResultsC6 | None

    def __init__(self, parent: SimulationPage) -> None:
        super().__init__(parent)

        self._xlim_cid: int | None = None

    def setupWidget(self) -> None:
        super().setupWidget()

        self._mappable = plt.cm.ScalarMappable(cmap=alphamagma, norm=Normalize(vmin=0, vmax=1))
        self._cbar = self.canvas.fig.colorbar(
            self._mappable,
            ax=self.canvas.ax,
            label=r"Admixture to perturbed state $\sum |V / \Delta E|$",
            aspect=60,
        )

    def plot(self, parameters: ParametersC6, results: ResultsC6) -> None:
        super().plot(parameters, results)

        if len(results.gaps_ghz) == 0:
            logger.warning("C6 Page: No contributions to plot.")
            return

        max_gap = np.max(np.abs(results.gaps_ghz))
        xlim = (-max_gap * 1.1, max_gap * 1.1) if max_gap > 0 else (-1, 1)

        self.draw_bars(parameters, results, xlim=xlim, ylim=None)

    def draw_bars(
        self,
        parameters: ParametersC6,
        results: ResultsC6,
        xlim: tuple[float, float],
        ylim: tuple[float, float] | None = None,
    ) -> None:
        ax = self.canvas.ax
        self.clear()
        self.parameters = parameters
        self.results = results

        self.bin_width_ghz = self.BIN_WIDTH_FRACTION * (xlim[1] - xlim[0])

        bin_to_stateinds: dict[int, list[int]] = {}
        bin_to_total: dict[int, float] = {}
        bin_to_admixture: dict[int, float] = {}
        for i, gap in enumerate(results.gaps_ghz):
            _bin = int(np.floor(gap / self.bin_width_ghz))
            contribution = results.contributions_c6[i]
            bin_to_stateinds.setdefault(_bin, []).append(i)
            bin_to_total[_bin] = bin_to_total.get(_bin, 0.0) + contribution
            bin_to_admixture[_bin] = bin_to_admixture.get(_bin, 0.0) + results.admixtures[i]

        show_status_tip(self, "Plotting contributions...")
        bins = sorted(bin_to_total)
        centers = [(b + 0.5) * self.bin_width_ghz for b in bins]
        heights = [bin_to_total[b] for b in bins]
        admixtures = np.array([bin_to_admixture[b] for b in bins])

        vmax = float(admixtures.max()) if len(admixtures) and admixtures.max() > 0 else 1.0
        self._mappable.set_norm(Normalize(vmin=0, vmax=vmax))
        self._cbar.update_normal(self._mappable)
        colors = self._mappable.to_rgba(admixtures)

        self.artist = ax.bar(centers, heights, width=self.bin_width_ghz * 0.9, color=colors)
        self._bar_data = {}
        for _bin, rect in zip(bins, self.artist.patches, strict=True):
            self._bar_data[_bin] = rect

        ax.scatter(centers, heights, color="black", s=10, zorder=5)

        ax.axhline(0, color="black", lw=0.5)
        ax.axvline(0, color="black", lw=0.5)
        ax.set_xlabel("Energy gap to state of interest (GHz)")
        ax.set_ylabel(r"Contribution to $C_6$ (GHz·μm⁶)")
        ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        else:
            ylim_max = np.max(np.abs(ax.get_ylim()))
            ax.set_ylim((-ylim_max, ylim_max))

        self.bin_to_stateinds = bin_to_stateinds
        self._xlim_cid = ax.callbacks.connect("xlim_changed", self._on_xlim_changed)

    def _on_xlim_changed(self, ax: mpl.axes.Axes) -> None:
        xlim = ax.get_xlim()

        new_bin_width_ghz = self.BIN_WIDTH_FRACTION * (xlim[1] - xlim[0])
        if abs(new_bin_width_ghz / self.bin_width_ghz - 1.0) < 0.1:
            return

        ylim = ax.get_ylim()
        parameters, results = self.parameters, self.results
        if parameters is None or results is None:
            return
        self.draw_bars(parameters, results, xlim=xlim, ylim=ylim)
        self.setup_annotations(parameters, results)

    def get_annotation_text(self, key: int) -> str:
        results = self.results
        if results is None:
            return "No data available."

        stateinds = self.bin_to_stateinds[key]
        stateinds = sorted(stateinds, key=lambda i: -abs(results.contributions_c6[i]))

        states = [results.kets[i] for i in stateinds]
        contributions = [results.contributions_c6[i] for i in stateinds]
        text = ""
        for state, contrib in zip(states, contributions, strict=True):
            text += f"  - {state}: {contrib:+.3f} GHz·μm⁶\n"

        return text

    def clear(self) -> None:
        super().clear()
        if self._xlim_cid is not None:
            self.canvas.ax.callbacks.disconnect(self._xlim_cid)
            self._xlim_cid = None


def fit_c3(x: NDArray[Any], /, e0: float, c3: float) -> NDArray[Any]:
    return e0 + c3 / x**3


def fit_c6(x: NDArray[Any], /, e0: float, c6: float) -> NDArray[Any]:
    return e0 + c6 / x**6


def fit_c3_c6(x: NDArray[Any], /, e0: float, c3: float, c6: float) -> NDArray[Any]:
    return e0 + c3 / x**3 + c6 / x**6
