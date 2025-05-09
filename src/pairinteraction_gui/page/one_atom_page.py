# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from typing import Any

from pairinteraction_gui.calculate.calculate_one_atom import ParametersOneAtom, ResultsOneAtom, calculate_one_atom
from pairinteraction_gui.config import (
    BasisConfigOneAtom,
    CalculationConfig,
    KetConfigOneAtom,
    SystemConfigOneAtom,
)
from pairinteraction_gui.page.base_page import CalculationPage

logger = logging.getLogger(__name__)


class OneAtomPage(CalculationPage):
    """Page for configuring and analyzing single atom systems."""

    title = "One Atom"
    tooltip = "Configure and analyze single atom systems"

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

    def update_plot(self, parameters: ParametersOneAtom, results: ResultsOneAtom) -> None:  # type: ignore[override]
        super().update_plot(parameters, results)
        self.add_short_labels(results)
        self.plotwidget.canvas.draw()

    def add_short_labels(
        self,
        results: ResultsOneAtom,
    ) -> None:
        if 0 not in results.state_labels:
            return

        ax = self.plotwidget.canvas.ax
        x_lim = ax.get_xlim()
        ax.set_xlim(x_lim[0] - (x_lim[1] - x_lim[0]) * 0.1, x_lim[1])

        used = set()
        for ket_label, energy in zip(results.state_labels[0], results.energies[0]):
            short_label = ket_label[1:-1]
            short_label = short_label.split(":", 1)[-1]
            components = short_label.split(",")
            short_label = ",".join(components[:-1])
            short_label = short_label.split("_", 1)[0]
            if short_label in used:
                continue
            used.add(short_label)
            self.plotwidget.canvas.ax.text(x_lim[0], energy, short_label, va="center", ha="right")

    def _get_export_notebook_template_name(self) -> str:
        return "one_atom.ipynb"

    def _get_export_replacements(self) -> dict[str, Any]:
        parameters = ParametersOneAtom.from_page(self)
        return parameters.to_replacement_dict()
