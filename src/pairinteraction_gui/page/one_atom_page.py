# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import pairinteraction._wrapped as pi
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
        self.system_atom_config = SystemConfigOneAtom(self)

    def calculate(self) -> None:
        super().calculate()

        self.fields = self.system_atom_config.get_fields()
        self.systems = self.system_atom_config.get_systems(0)
        self.ket = self.ket_config.get_ket_atom(0)
        kwargs: dict[str, Any] = {}
        if self.plotwidget.energy_range.isChecked():
            _energies = self.plotwidget.energy_range.values()
            kwargs["energy_range"] = (self.ket.get_energy("GHz") + v for v in _energies)
            kwargs["energy_unit"] = "GHz"
        if self.plotwidget.fast_mode.isChecked():
            kwargs["diagonalizer"] = "lapacke_evr"
            kwargs["float_type"] = "float32"

        pi.diagonalize(self.systems, **kwargs)

    def update_plot(self) -> None:
        ket_energy = self.ket.get_energy("GHz")
        energies = [system.get_eigenenergies("GHz") - ket_energy for system in self.systems]
        overlaps = [system.get_eigenbasis().get_overlaps(self.ket) for system in self.systems]

        x_values, xlabel = self.plotwidget._get_x_values_and_label_from_ranges(self.fields)

        self.plotwidget.plot(x_values, energies, overlaps, xlabel=xlabel)

        self.add_short_labels(energies)
        self.plotwidget.add_cursor(x_values[0], energies[0], self.systems[0])

        self.plotwidget.canvas.draw()

    def add_short_labels(
        self,
        energies: Sequence["NDArray[Any]"],
    ) -> None:
        ax = self.plotwidget.canvas.ax
        x_lim = ax.get_xlim()
        ax.set_xlim(x_lim[0] - (x_lim[1] - x_lim[0]) * 0.1, x_lim[1])

        basis0 = self.systems[0].get_eigenbasis()
        corresponding_kets = [basis0.get_corresponding_ket(i) for i in range(basis0.number_of_states)]

        used = set()
        l_dict: dict[float, str] = {0: "S", 1: "P", 2: "D", 3: "F", 4: "G", 5: "H"}
        for ket, energy in zip(corresponding_kets, energies[0]):
            if "mqdt" not in ket.species:
                short_label = f"{ket.n} {l_dict.get(ket.l, ket.l)}"
            else:
                short_label = f"{ket.n} L={ket.l:.1f}"
            if short_label in used:
                continue
            used.add(short_label)
            self.plotwidget.canvas.ax.text(x_lim[0], energy, short_label, va="center", ha="right")
