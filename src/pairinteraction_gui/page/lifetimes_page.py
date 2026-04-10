# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from pairinteraction_gui.calculate.calculate_lifetimes import ParametersLifetimes, calculate_lifetimes
from pairinteraction_gui.config import KetConfigLifetimes
from pairinteraction_gui.page.base_page import CalculationPage
from pairinteraction_gui.plotwidget.plotwidget import PlotLifetimes
from pairinteraction_gui.qobjects import show_status_tip

if TYPE_CHECKING:
    from collections.abc import Callable

    from pairinteraction_gui.calculate.calculate_lifetimes import ResultsLifetimes

logger = logging.getLogger(__name__)


class LifetimesPage(CalculationPage):
    """Page for calculating lifetimes."""

    title = "Lifetimes"
    tooltip = "Calculate the lifetimes and transition rates for a specified ket"

    ket_config: KetConfigLifetimes
    plotwidget: PlotLifetimes  # type: ignore[assignment]

    def setupWidget(self) -> None:
        super().setupWidget()
        show_status_tip(self, "Ready", timeout=1)

        # all attributes of instance BaseConfig will be added to the toolbox in postSetupWidget
        self.ket_config = KetConfigLifetimes(self)

    def _create_plot_widget(self) -> PlotLifetimes:  # type: ignore[override]
        return PlotLifetimes(self)

    def _get_export_actions(self) -> list[tuple[str, Callable[[], None]]]:
        return [("Export as PNG", self.export_png)]

    def calculate(self) -> tuple[ParametersLifetimes, ResultsLifetimes]:  # type: ignore[override]
        params = ParametersLifetimes(
            species=self.ket_config.get_species(),
            quantum_numbers=self.ket_config.get_quantum_numbers(),
            temperature=self.ket_config.get_temperature(),
        )
        results = calculate_lifetimes(params)
        return params, results

    def update_plot(self, parameters: ParametersLifetimes, results: ResultsLifetimes) -> None:  # type: ignore[override]
        super().update_plot(parameters, results)  # type: ignore[arg-type]
        self.ket_config.set_lifetime(results.lifetime)
        show_status_tip(self, "Finished updating plot. Tip: Click on a bar to see transition details.", logger=logger)
