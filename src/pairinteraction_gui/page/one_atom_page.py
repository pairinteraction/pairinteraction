# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

from pairinteraction_gui.calculate.calculate_one_atom import ParametersOneAtom, calculate_one_atom
from pairinteraction_gui.config import (
    BasisConfigOneAtom,
    KetConfigOneAtom,
    SystemConfigOneAtom,
)
from pairinteraction_gui.page.base_page import SimulationPage
from pairinteraction_gui.plotwidget.plotwidget import PlotEnergies

if TYPE_CHECKING:
    from numpy.typing import NDArray

logger = logging.getLogger(__name__)


class OneAtomPage(SimulationPage):
    """Page for configuring and analyzing single atom systems."""

    title = "One Atom"
    tooltip = "Configure and analyze single atom systems"

    plotwidget: PlotEnergies

    def setupWidget(self) -> None:
        self.plotwidget = PlotEnergies(self)
        self.layout().addWidget(self.plotwidget)
        super().setupWidget()

        # all attributes of instance BaseConfig will be added to the toolbox in postSetupWidget
        self.ket_config = KetConfigOneAtom(self)
        self.basis_config = BasisConfigOneAtom(self)
        self.system_config = SystemConfigOneAtom(self)

    def calculate(self) -> None:
        super().calculate()

        self.parameters = ParametersOneAtom.from_page(self)
        self.results = calculate_one_atom(self.parameters)

    def update_plot(self) -> None:
        energies = self.results.energies
        overlaps = self.results.ket_overlaps

        x_values = self.parameters.get_x_values()
        x_label = self.parameters.get_x_label()

        self.plotwidget.plot(x_values, energies, overlaps, x_label)

        self.add_short_labels(energies)
        self.plotwidget.add_cursor(x_values[0], energies[0], self.results.state_labels_0)

        self.plotwidget.canvas.draw()

    def add_short_labels(
        self,
        energies: Sequence["NDArray[Any]"],
    ) -> None:
        ax = self.plotwidget.canvas.ax
        x_lim = ax.get_xlim()
        ax.set_xlim(x_lim[0] - (x_lim[1] - x_lim[0]) * 0.1, x_lim[1])

        used = set()
        for ket_label, energy in zip(self.results.state_labels_0, energies[0]):
            short_label = ket_label[1:-1]
            short_label = short_label.split(":", 1)[-1]
            components = short_label.split(",")
            short_label = ",".join(components[:-1])
            short_label = short_label.split("_", 1)[0]
            if short_label in used:
                continue
            used.add(short_label)
            self.plotwidget.canvas.ax.text(x_lim[0], energy, short_label, va="center", ha="right")
