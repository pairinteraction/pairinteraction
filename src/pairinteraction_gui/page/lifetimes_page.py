# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import mplcursors
import numpy as np
from PySide6.QtGui import QPalette

from pairinteraction_gui.calculate.calculate_lifetimes import ParametersLifetimes, calculate_lifetimes
from pairinteraction_gui.config import KetConfigLifetimes
from pairinteraction_gui.page.base_page import CalculationPage
from pairinteraction_gui.plotwidget.plotwidget import PlotWidget
from pairinteraction_gui.qobjects import show_status_tip
from pairinteraction_gui.theme import theme_manager

if TYPE_CHECKING:
    from collections.abc import Callable

    from pairinteraction_gui.calculate.calculate_lifetimes import KetData, ResultsLifetimes

logger = logging.getLogger(__name__)


class LifetimesPage(CalculationPage):
    """Page for calculating lifetimes."""

    title = "Lifetimes"
    tooltip = "Calculate the lifetimes and transition rates for a specified ket"

    ket_config: KetConfigLifetimes
    plotwidget: PlotWidget  # type: ignore[assignment]

    def setupWidget(self) -> None:
        super().setupWidget()
        show_status_tip(self, "Ready", timeout=1)

        window_color = theme_manager.get_palette().color(QPalette.ColorRole.Window).name()

        self.plotwidget.canvas.fig.set_facecolor(window_color)
        self.plotwidget.canvas.fig.set_layout_engine(
            "constrained",
            w_pad=0.2,
            h_pad=0.2,
            wspace=0.0,
            hspace=0.0,
        )

        # all attributes of instance BaseConfig will be added to the toolbox in postSetupWidget
        self.ket_config = KetConfigLifetimes(self)

    def _create_plot_widget(self) -> PlotWidget:  # type: ignore[override]
        return PlotWidget(self)

    def _get_export_actions(self) -> list[tuple[str, Callable[[], None]]]:
        return [("Export as PNG", self.export_png)]

    def calculate(self) -> tuple[ParametersLifetimes, ResultsLifetimes]:  # type: ignore[override]
        params = ParametersLifetimes(
            species=self.ket_config.get_species(),
            quantum_numbers=self.ket_config.get_quantum_numbers(),
            temperature=self.ket_config.get_temperature(),
        )
        results = calculate_lifetimes(params)
        return params, results

    def update_plot(self, parameters: ParametersLifetimes, results: ResultsLifetimes) -> None:  # type: ignore[override]
        ax = self.plotwidget.canvas.ax
        ax.clear()

        n_list = np.arange(0, np.max([s.n for s in results.kets_bbr + results.kets_sp] + [0]) + 1)
        sorted_rates: dict[str, dict[int, list[tuple[KetData, float]]]] = {}
        for key, kets, rates in [
            ("BBR", results.kets_bbr, results.transition_rates_bbr),
            ("SP", results.kets_sp, results.transition_rates_sp),
        ]:
            sorted_rates[key] = {n: [] for n in n_list}
            for i, s in enumerate(kets):
                show_status_tip(self, "Preparing transition rates...")
                sorted_rates[key][s.n].append((s, rates[i]))
        self.sorted_rates = sorted_rates

        show_status_tip(self, "Plotting transition rates...")
        rates_summed = {key: [sum(rates for _, rates in sorted_rates[key][n]) for n in n_list] for key in sorted_rates}
        bar_sp = ax.bar(n_list, rates_summed["SP"], label="Spontaneous Decay", color="blue", alpha=0.8)
        bar_bbr = ax.bar(n_list, rates_summed["BBR"], label="Black Body Radiation", color="red", alpha=0.8)
        self.artists = (bar_sp, bar_bbr)
        ax.legend()

        ax.set_xlabel("Principal Quantum Number $n$")
        ax.set_ylabel(r"Transition Rates (1 / ms)")

        show_status_tip(self, "Adding transition rate annotations...")
        self.add_cursor()

        self.plotwidget.canvas.draw()
        self.plotwidget.navigation_toolbar.reset_home_view()
        self.ket_config.set_lifetime(results.lifetime)

    def add_cursor(self) -> None:
        """Add interactive cursor to the plot."""
        # Remove any existing cursors to avoid duplicates
        if hasattr(self, "mpl_cursor"):
            if hasattr(self.mpl_cursor, "remove"):  # type: ignore
                self.mpl_cursor.remove()  # type: ignore
            del self.mpl_cursor  # type: ignore

        self.mpl_cursor = mplcursors.cursor(
            self.artists,
            hover=mplcursors.HoverMode.Transient,
            annotation_kwargs={
                "bbox": {"boxstyle": "round,pad=0.5", "fc": "white", "alpha": 0.9, "ec": "gray"},
                "arrowprops": {"arrowstyle": "->", "connectionstyle": "arc3", "color": "gray"},
            },
        )

        @self.mpl_cursor.connect("add")
        def on_add(sel: mplcursors.Selection) -> None:
            label = sel.artist.get_label()
            x, y, width, height = sel.artist[sel.index].get_bbox().bounds

            n = round(x + width / 2)
            key = "BBR" if "Black Body" in label else "SP"
            state_text = "\n".join(f"{s}: {r:.5f}/ms" for (s, r) in self.sorted_rates[key][n])
            text = f"{label} to n={n}:\n{state_text}"

            sel.annotation.set(text=text, position=(0, 20), anncoords="offset points")
            sel.annotation.xy = (x + width / 2, y + height / 2)
