# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from pairinteraction_gui.calculate.calculate_one_atom import ParametersOneAtom, calculate_one_atom
from pairinteraction_gui.config import (
    BasisConfigOneAtom,
    CalculationConfig,
    KetConfigOneAtom,
    SystemConfigOneAtom,
)
from pairinteraction_gui.page.base_page import CalculationPage

if TYPE_CHECKING:
    from pairinteraction_gui.calculate.calculate_one_atom import ResultsOneAtom

logger = logging.getLogger(__name__)


class OneAtomPage(CalculationPage):
    """Page for configuring and analyzing single-atom systems."""

    title = "One\nAtom"
    tooltip = "Configure and analyze single-atom systems"

    def setupWidget(self) -> None:
        super().setupWidget()

        # all attributes of instance BaseConfig will be added to the toolbox in postSetupWidget
        self.ket_config = KetConfigOneAtom(self)
        self.basis_config = BasisConfigOneAtom(self)
        self.system_config = SystemConfigOneAtom(self)
        self.calculation_config = CalculationConfig(self)

        self.ket_config.signal_species_changed.connect(self.basis_config.on_species_changed)
        self.ket_config.signal_species_changed.connect(self.plotwidget.clear)

    def calculate(self) -> tuple[ParametersOneAtom, ResultsOneAtom]:
        parameters = ParametersOneAtom.from_page(self)
        results = calculate_one_atom(parameters)
        return parameters, results

    def _get_export_notebook_template_name(self) -> str:
        return "one_atom.ipynb"

    def _get_export_replacements(self) -> dict[str, str]:
        parameters = ParametersOneAtom.from_page(self)
        return parameters.to_replacement_dict()
