import logging
from typing import Any, Union

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

    _export_notebook_template = "two_atoms.ipynb"

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
                    energy=(ket_pair_energy - delta_energy,
                     + delta_energy),
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
        energies = [system.get_eigenvalues("GHz") - ket_pair_energy for system in self.systems]
        overlaps = [system.get_eigenbasis().get_overlaps(self.kets) for system in self.systems]

        ranges = {**self.fields, "distance": self.distance, "angle": self.angle}
        x_values, xlabel = self.plotwidget._get_x_values_and_label_from_ranges(ranges)

        self.plotwidget.plot(x_values, energies, overlaps, xlabel=xlabel)

        self.plotwidget.add_cursor(x_values[-1], energies[-1], self.systems[-1])

        self.plotwidget.canvas.draw()

    def _get_export_replacements(self) -> dict[str, Any]:
        if not all(f.is_constant() for f in self.fields.values()):
            raise NotImplementedError("Exporting is currently only supported for constant fields.")

        replacements: dict[str, Union[str, float]] = {}
        for i in range(2):
            ket = self.ket_config.get_ket_atom(i)
            qns = self.ket_config.get_quantum_numbers(i)
            restrict_qns = self.basis_config.get_restrict_qns(i)

            replacements.update({
                f"$SPECIES{i+1}": f'"{ket.species}"',
                f"$QUANTUM_NUMBERS_RESTRICTIONS{i+1}": ", ".join(f"{k}={v}" for k, v in restrict_qns.items()),
                f"$QUANTUM_NUMBERS{i+1}": ", ".join(f"{k}={v}" for k, v in qns.items()),
            })

        ranges = {**self.fields, "distance": self.distance, "angle": self.angle}
        replacements.update({
            "$EX_VALUE": ranges["Ex"].min_value,
            "$EY_VALUE": ranges["Ey"].min_value,
            "$EZ_VALUE": ranges["Ez"].min_value,
            "$BX_VALUE": ranges["Bx"].min_value,
            "$BY_VALUE": ranges["By"].min_value,
            "$BZ_VALUE": ranges["Bz"].min_value,
        })
        replacements.update({
            "$STEPS": self.system_pair_config.steps_spinbox.value(),
            "$DISTANCE_MIN": ranges["distance"].min_value,
            "$DISTANCE_MAX": ranges["distance"].max_value,
            "$ANGLE_MIN": ranges["angle"].min_value,
            "$ANGLE_MAX": ranges["angle"].max_value,
        })
        replacements["$MULTIPOLE_ORDER"] = self.system_pair_config.get_order()

        # PAIR_ENERGY_DELTA
        replacements["$PAIR_ENERGY_DELTA"] = self.basis_config.delta_pair_energy.value()

        # Add diagonalization kwargs if needed
        diag_kwargs: dict[str, Any] = {}
        # if self.plotwidget.energy_range.isChecked():
        #     _energies = self.plotwidget.energy_range.values()
        #     diag_kwargs["energy_range"] = tuple(ket.get_energy("GHz") + v for v in _energies)
        #     diag_kwargs["energy_unit"] = "GHz"
        if self.plotwidget.fast_mode.isChecked():
            diag_kwargs["diagonalizer"] = "lapacke_evr"
            diag_kwargs["float_type"] = "float32"
        replacements["$DIAGONALIZE_KWARGS"] = ""
        if diag_kwargs:
            replacements["$DIAGONALIZE_KWARGS"] = ", " + ", ".join(f"{k}={v!r}" for k, v in diag_kwargs.items())


        # Handle x-axis values and label
        _, xlabel = self.plotwidget._get_x_values_and_label_from_ranges(ranges)
        replacements["$X_VALUES"] = xlabel.split(" ")[0]
        replacements["$XLABEL"] = f'r"{xlabel}"'

        isreal = ranges["By"].is_zero() and ranges["Ey"].is_zero()
        replacements["$PI_DTYPE"] = "real" if isreal else "complex"

        return replacements
