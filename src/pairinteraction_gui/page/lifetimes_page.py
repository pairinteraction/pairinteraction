# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from typing import TYPE_CHECKING

import mplcursors
import numpy as np

from pairinteraction_gui.config import KetConfigLifetimes
from pairinteraction_gui.page.base_page import SimulationPage
from pairinteraction_gui.plotwidget.plotwidget import PlotWidget
from pairinteraction_gui.qobjects import show_status_tip

if TYPE_CHECKING:
    import pairinteraction.real as pi

logger = logging.getLogger(__name__)


class LifetimesPage(SimulationPage):
    """Page for calculating lifetimes."""

    title = "Lifetimes"
    tooltip = "Calculate the lifetimes and transition rates for a specified ket."

    ket_config: KetConfigLifetimes

    def setupWidget(self) -> None:
        self.plotwidget = PlotWidget(self)
        self.layout().addWidget(self.plotwidget)
        super().setupWidget()

        show_status_tip(self, "Ready", timeout=1)

        # all attributes of instance BaseConfig will be added to the toolbox in postSetupWidget
        self.ket_config = KetConfigLifetimes(self)

    def calculate(self) -> None:
        # since this calculation is rather fast, we can just call it directly and dont have to use a different process
        ket = self.ket_config.get_ket_atom(0)
        temperature = self.ket_config.get_temperature()
        self.kets_sp, self.transition_rates_sp = ket.get_spontaneous_transition_rates(unit="1/ms")
        self.kets_bbr, self.transition_rates_bbr = ket.get_black_body_transition_rates(temperature, "K", unit="1/ms")

    def update_plot(self) -> None:
        ax = self.plotwidget.canvas.ax
        ax.clear()

        n_list = np.arange(0, np.max([s.n for s in self.kets_bbr]) + 1)
        sorted_rates: dict[str, dict[int, list[tuple[pi.KetAtom, float]]]] = {}
        for key, kets, rates in [
            ("BBR", self.kets_bbr, self.transition_rates_bbr),
            ("SP", self.kets_sp, self.transition_rates_sp),
        ]:
            sorted_rates[key] = {n: [] for n in n_list}
            for i, s in enumerate(kets):
                sorted_rates[key][s.n].append((s, rates[i]))
        self.sorted_rates = sorted_rates

        rates_summed = {key: [sum(rates for _, rates in sorted_rates[key][n]) for n in n_list] for key in sorted_rates}
        bar_sp = ax.bar(n_list, rates_summed["SP"], label="Spontaneous Decay", color="blue", alpha=0.8)
        bar_bbr = ax.bar(n_list, rates_summed["BBR"], label="Black Body Radiation", color="red", alpha=0.8)
        self.artists = (bar_sp, bar_bbr)
        ax.legend()

        ax.set_xlabel("Principal Quantum Number $n$")
        ax.set_ylabel(r"Transition Rates (1 / ms)")

        self.add_cursor()

        self.plotwidget.canvas.draw()

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
