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

    def calculate(self) -> tuple[ParametersLifetimes, ResultsLifetimes]:  # type: ignore[override]
        params = ParametersLifetimes.from_page(self)
        results = calculate_lifetimes(params)
        return params, results

    def update_plot(self, parameters: ParametersLifetimes, results: ResultsLifetimes) -> None:  # type: ignore[override]
        super().update_plot(parameters, results)  # type: ignore[arg-type]
        self.ket_config.set_lifetime(results.lifetime)
        show_status_tip(self, "Finished updating plot. Tip: Click on a bar to see transition details.", logger=logger)

    def _get_export_notebook_template_name(self) -> str:
        return "lifetimes.ipynb"

    def _get_export_replacements(self) -> dict[str, str]:
        parameters = ParametersLifetimes.from_page(self)
        return parameters.to_replacement_dict()
