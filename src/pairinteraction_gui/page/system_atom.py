import logging
from typing import Any

import pairinteraction._wrapped as pi
from pairinteraction_gui.config import BasisAtomConfig, KetAtomConfig, SystemAtomConfig
from pairinteraction_gui.page.base_page import SimulationPage
from pairinteraction_gui.plotwidget.plotwidget import PlotEnergies

logger = logging.getLogger(__name__)


class SystemAtomPage(SimulationPage):
    """Page for configuring and analyzing single atom systems."""

    title = "One Atom"
    tooltip = "Configure and analyze single atom systems"

    plotwidget: PlotEnergies

    def setupWidget(self) -> None:
        self.plotwidget = PlotEnergies(self)
        self.layout().addWidget(self.plotwidget)
        super().setupWidget()

        # all attributes of instance BaseConfig will be added to the toolbox in postSetupWidget
        self.ket_config = KetAtomConfig(self)
        self.basis_config = BasisAtomConfig(self)
        self.system_atom_config = SystemAtomConfig(self)

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
        energies = [system.get_eigenvalues("GHz") - ket_energy for system in self.systems]
        overlaps = [system.get_eigenbasis().get_overlaps(self.ket) for system in self.systems]

        self.plotwidget.plot(energies, self.fields, overlaps)

        self.plotwidget.canvas.draw()
