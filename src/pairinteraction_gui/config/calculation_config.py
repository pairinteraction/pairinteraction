# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later


from PySide6.QtWidgets import QLabel

from pairinteraction_gui.config.base_config import BaseConfig
from pairinteraction_gui.qobjects import Item, QnItemInt, RangeItem


class CalculationConfig(BaseConfig):
    """Section for configuring the system parameters."""

    margin = (5, 15, 5, 5)
    spacing = 10

    title = "Calculation Options"

    def setupWidget(self) -> None:
        self.layout().addWidget(QLabel("<b>Calculation steps</b>"))
        self.steps = QnItemInt(
            self,
            "Number of steps",
            unit="",
            vmin=1,
            vmax=9999,
            vdefault=100,
            tooltip="Number of steps for the calculation",
        )
        self.layout().addWidget(self.steps)

        self.layout().addWidget(QLabel("<b>Fast mode</b>"))
        self.fast_mode = Item(self, "Use fast calculation mode", checked=True)
        self.layout().addWidget(self.fast_mode)

        self.layout().addWidget(QLabel("<b>State Annotations</b> (click on the diamond shaped markers)"))
        self.number_state_labels = QnItemInt(
            self,
            "Annotate",
            unit="state labels",
            vdefault=10,
            tooltip="Number of steps, for which the overlap is calculated and the corresponding state labels"
            " can be shown by clicking on the info circles.",
        )
        self.layout().addWidget(self.number_state_labels)

        self.layout().addWidget(QLabel("<b>Energy Range</b>"))
        self.layout().addWidget(QLabel("Calculate eigenenergies in the range"))
        self.energy_range = RangeItem(
            self,
            "from",
            vdefaults=(-80, 80),
            unit="GHz",
            checked=False,
            tooltip_label="energy",
        )
        self.layout().addWidget(self.energy_range)
