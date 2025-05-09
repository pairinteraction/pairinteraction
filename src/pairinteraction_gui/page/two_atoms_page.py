# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from typing import Any

from pairinteraction_gui.calculate.calculate_two_atoms import ParametersTwoAtoms, ResultsTwoAtoms, calculate_two_atoms
from pairinteraction_gui.config import (
    BasisConfigTwoAtoms,
    CalculationConfig,
    KetConfigTwoAtoms,
    SystemConfigTwoAtoms,
)
from pairinteraction_gui.page.base_page import CalculationPage

logger = logging.getLogger(__name__)


class TwoAtomsPage(CalculationPage):
    """Page for configuring and analyzing pair systems."""

    title = "Two Atoms"
    tooltip = "Configure and analyze pair systems"

    def setupWidget(self) -> None:
        super().setupWidget()

        # all attributes of instance BaseConfig will be added to the toolbox in postSetupWidget
        self.ket_config = KetConfigTwoAtoms(self)
        self.basis_config = BasisConfigTwoAtoms(self)
        self.system_config = SystemConfigTwoAtoms(self)
        self.calculation_config = CalculationConfig(self)

        # Set some better defaults for the two atoms page
        self.calculation_config.number_state_labels.setValue(5)
        self.calculation_config.number_state_labels.setChecked(False)
        self.calculation_config.energy_range.setValues(-0.5, 0.5)
        self.calculation_config.energy_range.setChecked(False)

        self.ket_config.signal_species_changed.connect(self.basis_config.on_species_changed)
        self.ket_config.signal_species_changed.connect(self.plotwidget.clear)

    def before_calculate(self) -> None:
        self.basis_config.clear_basis_pair_label()
        return super().before_calculate()

    def calculate(self) -> tuple[ParametersTwoAtoms, ResultsTwoAtoms]:
        parameters = ParametersTwoAtoms.from_page(self)
        results = calculate_two_atoms(parameters)
        return parameters, results

    def update_plot(self, parameters: ParametersTwoAtoms, results: ResultsTwoAtoms) -> None:  # type: ignore[override]
        super().update_plot(parameters, results)

        if results.basis_0_label is not None:
            self.basis_config.update_basis_pair_label(results.basis_0_label)

    def _get_export_replacements(self) -> dict[str, Any]:
        parameters = ParametersTwoAtoms.from_page(self)
        return parameters.to_replacement_dict()

    def _get_export_notebook_template_name(self) -> str:
        ranges = self.system_config.get_ranges_dict()
        if all(v[0] == v[-1] for k, v in ranges.items() if k in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]):
            return "two_atoms.ipynb"
        return "two_atoms_variable_fields.ipynb"
