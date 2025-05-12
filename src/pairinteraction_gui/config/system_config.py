# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Literal, Union

import numpy as np
from PySide6.QtWidgets import (
    QLabel,
)

from pairinteraction_gui.config.base_config import BaseConfig
from pairinteraction_gui.qobjects.item import QnItemInt, RangeItem

if TYPE_CHECKING:
    from pairinteraction_gui.page import OneAtomPage, TwoAtomsPage


RangesKeys = Literal["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Distance", "Angle"]


class SystemConfig(BaseConfig):
    """Section for configuring the system parameters."""

    margin = (5, 15, 5, 5)
    spacing = 10

    title = "System"
    page: Union["OneAtomPage", "TwoAtomsPage"]

    def setupEField(self) -> None:
        efield_label = QLabel("<b>Electric field</b>")
        self.layout().addWidget(efield_label)

        self.Ex = RangeItem(self, "Ex", unit="V/cm", tooltip_label="electric field in x-direction")
        self.Ey = RangeItem(self, "Ey", unit="V/cm", tooltip_label="electric field in y-direction")
        self.Ez = RangeItem(self, "Ez", unit="V/cm", tooltip_label="electric field in z-direction")

        self.layout().addWidget(self.Ex)
        self.layout().addWidget(self.Ey)
        self.layout().addWidget(self.Ez)

    def setupBField(self) -> None:
        bfield_label = QLabel("<b>Magnetic field</b>")
        self.layout().addWidget(bfield_label)

        self.Bx = RangeItem(self, "Bx", unit="Gauss", tooltip_label="magnetic field in x-direction")
        self.By = RangeItem(self, "By", unit="Gauss", tooltip_label="magnetic field in y-direction")
        self.Bz = RangeItem(self, "Bz", unit="Gauss", tooltip_label="magnetic field in z-direction")

        self.layout().addWidget(self.Bx)
        self.layout().addWidget(self.By)
        self.layout().addWidget(self.Bz)

    def get_ranges_dict(self) -> dict[RangesKeys, list[float]]:
        """Return the electric and magnetic field ranges."""
        steps = self.page.calculation_config.steps.value()
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


class SystemConfigTwoAtoms(SystemConfig):
    page: "TwoAtomsPage"

    def setupWidget(self) -> None:
        self.setupEField()
        self.setupBField()
        self.setupDistance()
        self.setupAngle()
        self.setupOrder()

    def setupDistance(self) -> None:
        label = QLabel("<b>Distance</b>")
        self.layout().addWidget(label)

        self.distance = RangeItem(self, "Distance", vdefaults=(3, 8), vrange=(0, np.inf), unit="<span>&mu;m</span>")
        self.layout().addWidget(self.distance)

    def setupAngle(self) -> None:
        label = QLabel("<b>Angle</b> (0° = z-axis, 90° = x-axis)")
        self.layout().addWidget(label)

        self.angle = RangeItem(self, "Angle", vdefaults=(0, 0), vrange=(0, 360), unit="degree")
        self.layout().addWidget(self.angle)

    def setupOrder(self) -> None:
        self.layout().addWidget(QLabel("<b>Multipole expansion order</b>"))
        self.order = QnItemInt(
            self,
            "Multipole order",
            vmin=3,
            vmax=5,
            vdefault=3,
            tooltip="Select the order of the multipole expansion",
            checkable=False,
        )
        self.layout().addWidget(self.order)

    def _get_all_ranges(self) -> list[RangeItem]:
        """Return all range items."""
        return [*super()._get_all_ranges(), self.distance, self.angle]
