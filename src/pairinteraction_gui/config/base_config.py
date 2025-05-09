# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING

from pairinteraction_gui.qobjects import WidgetV

if TYPE_CHECKING:
    from pairinteraction_gui.page import SimulationPage


class BaseConfig(WidgetV):
    """Base class for configurations."""

    title: str

    def __init__(self, parent: "SimulationPage") -> None:
        """Initialize the base section."""
        self.page = parent
        super().__init__(parent)

    def postSetupWidget(self) -> None:
        self.layout().addStretch(1)
