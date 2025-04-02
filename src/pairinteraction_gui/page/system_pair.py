# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from typing import Any

import numpy as np

from pairinteraction import (
    complex as pi_complex,
    real as pi_real,
)
from pairinteraction_gui.config import BasisPairConfig, KetPairConfig, SystemPairConfig
from pairinteraction_gui.page.base_page import SimulationPage
from pairinteraction_gui.plotwidget.plotwidget import PlotEnergies

logger = logging.getLogger(__name__)


class SystemPairPage(SimulationPage):
    """Page for configuring and analyzing pair systems."""

    title = "Two Atoms"
    tooltip = "Configure and analyze pair systems"

    plotwidget: PlotEnergies

    def setupWidget(self) -> None:
        self.plotwidget = PlotEnergies(self)
        self.layout().addWidget(self.plotwidget)
        super().setupWidget()

        # all attributes of instance BaseConfig will be added to the toolbox in postSetupWidget
        self.ket_config = KetPairConfig(self)
        self.basis_config = BasisPairConfig(self)
        self.system_pair_config = SystemPairConfig(self)

    def calculate(self) -> None:
        super().calculate()
        self.basis_config.clear_basis_pair_label()

        self.fields = self.system_pair_config.get_fields()
        self.distance = self.system_pair_config.get_distance()
        self.angle = self.system_pair_config.get_angle()
        self.kets = [self.ket_config.get_ket_atom(i) for i in range(2)]
        delta_energy = self.basis_config.delta_pair_energy.value()

        kwargs: dict[str, Any] = {}
        if self.plotwidget.fast_mode.isChecked():
            kwargs["diagonalizer"] = "lapacke_evr"
            kwargs["float_type"] = "float32"

        isreal = self.fields["Ey"].is_zero() and self.fields["By"].is_zero()
        pi = pi_real if isreal else pi_complex

        steps = self.distance.steps
        if all(f.is_constant() for f in self.fields.values()):
            system_atoms = [self.system_pair_config.get_systems(i)[0] for i in range(2)]
            pi.diagonalize(system_atoms, **kwargs)
            ket_pair_energy = sum(
                system.get_corresponding_energy(ket, "GHz") for ket, system in zip(self.kets, system_atoms)
            )
            basis_pair = pi.BasisPair(
                system_atoms, energy=(ket_pair_energy - delta_energy, ket_pair_energy + delta_energy), energy_unit="GHz"
            )
            basis_pair_list = [basis_pair]
        else:
            system_atoms_list = [self.system_pair_config.get_systems(i) for i in range(2)]
            pi.diagonalize(system_atoms_list[0] + system_atoms_list[1], **kwargs)
            basis_pair_list = []
            for step in range(steps):
                system_atoms = [system_atoms_list[0][step], system_atoms_list[1][step]]
                ket_pair_energy = sum(
                    system.get_corresponding_energy(ket, "GHz") for ket, system in zip(self.kets, system_atoms)
                )
                basis_pair = pi.BasisPair(
                    system_atoms,
                    energy=(ket_pair_energy - delta_energy, ket_pair_energy + delta_energy),
                    energy_unit="GHz",
                )
                basis_pair_list.append(basis_pair)

        self.systems = []
        for step in range(steps):
            basis_pair = basis_pair_list[step % len(basis_pair_list)]
            system = pi.SystemPair(basis_pair)
            system.set_order(self.system_pair_config.get_order())
            if not np.isinf(self.distance[step]):
                system.set_distance(self.distance[step], self.angle[step], unit="micrometer")
            self.systems.append(system)

        if self.plotwidget.energy_range.isChecked():
            ket_pair_energy = sum(ket.get_energy("GHz") for ket in self.kets)
            _energies = self.plotwidget.energy_range.values()
            kwargs["energy_range"] = (ket_pair_energy + v for v in _energies)
            kwargs["energy_unit"] = "GHz"

        self.basis_config.update_basis_pair_label(basis_pair_list[0])
        pi.diagonalize(self.systems, **kwargs)

    def update_plot(self) -> None:
        ket_pair_energy = sum(ket.get_energy("GHz") for ket in self.kets)
        energies = [system.get_eigenenergies("GHz") - ket_pair_energy for system in self.systems]
        overlaps = [system.get_eigenbasis().get_overlaps(self.kets) for system in self.systems]

        ranges = {**self.fields, "distance": self.distance, "angle": self.angle}
        x_values, xlabel = self.plotwidget._get_x_values_and_label_from_ranges(ranges)

        self.plotwidget.plot(x_values, energies, overlaps, xlabel=xlabel)

        self.plotwidget.add_cursor(x_values[-1], energies[-1], self.systems[-1])

        self.plotwidget.canvas.draw()
