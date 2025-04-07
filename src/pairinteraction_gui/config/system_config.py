# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Literal, Union

import numpy as np
from PySide6.QtWidgets import (
    QLabel,
)

from pairinteraction_gui.config.base_config import BaseConfig
from pairinteraction_gui.qobjects import IntSpinBox, WidgetForm
from pairinteraction_gui.qobjects.item import RangeItem

if TYPE_CHECKING:
    from pairinteraction_gui.page import OneAtomPage, TwoAtomsPage


RangesKeys = Literal["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Distance", "Angle"]


class SystemConfig(BaseConfig):
    """Section for configuring the system parameters."""

    margin = (5, 15, 5, 5)
    spacing = 10

    title = "System"
    page: Union["OneAtomPage", "TwoAtomsPage"]

    def postSetupWidget(self) -> None:
        self.layout().addStretch()

    def setupEField(self) -> None:
        efield_label = QLabel("<b>Electric field</b>")
        self.layout().addWidget(efield_label)

        self.Ex = RangeItem(self, "Ex", (-999, 999), unit="V/cm", tooltip_label="electric field in x-direction")
        self.Ey = RangeItem(self, "Ey", (-999, 999), unit="V/cm", tooltip_label="electric field in y-direction")
        self.Ez = RangeItem(self, "Ez", (-999, 999), unit="V/cm", tooltip_label="electric field in z-direction")

        self.layout().addWidget(self.Ex)
        self.layout().addWidget(self.Ey)
        self.layout().addWidget(self.Ez)

    def setupBField(self) -> None:
        bfield_label = QLabel("<b>Magnetic field</b>")
        self.layout().addWidget(bfield_label)

        self.Bx = RangeItem(self, "Bx", (-999, 999), unit="Gauss", tooltip_label="magnetic field in x-direction")
        self.By = RangeItem(self, "By", (-999, 999), unit="Gauss", tooltip_label="magnetic field in y-direction")
        self.Bz = RangeItem(self, "Bz", (-999, 999), unit="Gauss", tooltip_label="magnetic field in z-direction")

        self.layout().addWidget(self.Bx)
        self.layout().addWidget(self.By)
        self.layout().addWidget(self.Bz)

    def setupSteps(self) -> None:
        steps_label = QLabel("<b>Calculation steps</b>")
        self.layout().addWidget(steps_label)

        steps_widget = WidgetForm(margin=(0, 0, 0, 0))
        self.steps_spinbox = IntSpinBox(vmin=1, vmax=10000, vdefault=100, tooltip="Number of steps for the calculation")
        steps_widget.layout().addRow("Number of steps", self.steps_spinbox)
        self.layout().addWidget(steps_widget)

    def get_ranges_dict(self) -> dict[RangesKeys, list[float]]:
        """Return the electric and magnetic field ranges."""
        steps = self.steps_spinbox.value()
        all_ranges = self._get_all_ranges()
        ranges_min_max: dict[str, tuple[float, float]] = {
            item.label.text(): item.values() for item in all_ranges if item.isChecked()
        }
        if len(ranges_min_max) == 0:
            ranges_min_max["Bz"] = (0, 0)

        return {key: np.linspace(value[0], value[1], steps).tolist() for key, value in ranges_min_max.items()}  # type: ignore [misc]

    def _get_all_ranges(self) -> list[RangeItem]:
        """Return all range items."""
        return [self.Ex, self.Ey, self.Ez, self.Bx, self.By, self.Bz]


class SystemConfigOneAtom(SystemConfig):
    page: "OneAtomPage"

    def setupWidget(self) -> None:
        self.setupEField()
        self.setupBField()
        self.setupSteps()


class SystemConfigTwoAtoms(SystemConfig):
    page: "TwoAtomsPage"

    def setupWidget(self) -> None:
        self.setupEField()
        self.setupBField()
        self.setupDistance()
        self.setupAngle()
        self.setupOrder()
        self.setupSteps()

    def setupDistance(self) -> None:
        label = QLabel("<b>Distance</b>")
        self.layout().addWidget(label)

        self.distance = RangeItem(self, "Distance", (0, 999), (3, 8), unit="<span>&mu;m</span>")
        self.layout().addWidget(self.distance)

    def setupAngle(self) -> None:
        label = QLabel("<b>Angle</b> (0° = z-axis, 90° = x-axis)")
        self.layout().addWidget(label)

        self.angle = RangeItem(self, "Angle", (0, 360), (0, 0), unit="degree")
        self.layout().addWidget(self.angle)

    def setupOrder(self) -> None:
        steps_label = QLabel("<b>Multipole expansion order</b>")
        self.layout().addWidget(steps_label)

        steps_widget = WidgetForm(margin=(0, 0, 0, 0))
        self.order = IntSpinBox(self, vmin=3, vmax=5, vdefault=3, tooltip="Select the order of the multipole expansion")
        steps_widget.layout().addRow("Multipole order", self.order)
        self.layout().addWidget(steps_widget)

    def _get_all_ranges(self) -> list[RangeItem]:
        """Return all range items."""
        return [*super()._get_all_ranges(), self.distance, self.angle]
