# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from pairinteraction_gui.calculate.calculate_c6 import ParametersC6, calculate_c6
from pairinteraction_gui.config import BasisConfigC6, KetConfigC6, SystemConfigC6
from pairinteraction_gui.page.base_page import CalculationPage
from pairinteraction_gui.plotwidget.plotwidget import PlotC6
from pairinteraction_gui.qobjects import show_status_tip

if TYPE_CHECKING:
    from collections.abc import Callable

    from pairinteraction_gui.calculate.calculate_c6 import ResultsC6

logger = logging.getLogger(__name__)


class C6Page(CalculationPage):
    """Page for calculating the C6 coefficient."""

    title = "C6"
    tooltip = "Calculate the C6 coefficient for a specified pair of kets"

    ket_config: KetConfigC6
    basis_config: BasisConfigC6
    system_config: SystemConfigC6
    plotwidget: PlotC6  # type: ignore[assignment]

    def setupWidget(self) -> None:
        super().setupWidget()
        show_status_tip(self, "Ready", timeout=1)

        # all attributes of instance BaseConfig will be added to the toolbox in postSetupWidget
        self.ket_config = KetConfigC6(self)
        self.basis_config = BasisConfigC6(self)
        self.system_config = SystemConfigC6(self)
        self.calculation_config = None

        self.ket_config.signal_species_changed.connect(self.basis_config.on_species_changed)
        self.ket_config.signal_species_changed.connect(self.plotwidget.clear)

    def _create_plot_widget(self) -> PlotC6:  # type: ignore[override]
        return PlotC6(self)

    def before_calculate(self) -> None:
        self.basis_config.clear_basis_pair_label()
        return super().before_calculate()

    def calculate(self) -> tuple[ParametersC6, ResultsC6]:  # type: ignore[override]
        params = ParametersC6.from_page(self)
        results = calculate_c6(params)
        return params, results

    def update_plot(self, parameters: ParametersC6, results: ResultsC6) -> None:  # type: ignore[override]
        super().update_plot(parameters, results)  # type: ignore[arg-type]
        self.ket_config.set_c6(results.c6)
        show_status_tip(self, "Finished updating plot. Tip: Click on a bar to see contributing kets.", logger=logger)

    def _get_export_actions(self) -> list[tuple[str, Callable[[], None]]]:
        return [("Export as PNG", self.export_png)]
