# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging

from pairinteraction_gui.calculate.calculate_two_atoms import ParametersTwoAtoms, ResultsTwoAtoms, calculate_two_atoms
from pairinteraction_gui.config import (
    BasisConfigTwoAtoms,
    KetConfigTwoAtoms,
    SystemConfigTwoAtoms,
)
from pairinteraction_gui.page.base_page import CalculationPage
from pairinteraction_gui.plotwidget.plotwidget import PlotEnergies

logger = logging.getLogger(__name__)


class TwoAtomsPage(CalculationPage):
    """Page for configuring and analyzing pair systems."""

    title = "Two Atoms"
    tooltip = "Configure and analyze pair systems"

    def setupWidget(self) -> None:
        self.plotwidget = PlotEnergies(self)
        self.layout().addWidget(self.plotwidget)
        super().setupWidget()

        # all attributes of instance BaseConfig will be added to the toolbox in postSetupWidget
        self.ket_config = KetConfigTwoAtoms(self)
        self.basis_config = BasisConfigTwoAtoms(self)
        self.system_config = SystemConfigTwoAtoms(self)

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
