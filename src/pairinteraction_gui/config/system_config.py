# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
from PySide6.QtWidgets import (
    QLabel,
)

from pairinteraction_gui.config.base_config import BaseConfig
from pairinteraction_gui.qobjects import CheckboxItem, ParameterItem, ParameterItemRange, QnItemInt

if TYPE_CHECKING:
    from pairinteraction_gui.page import C6Page, OneAtomPage, TwoAtomsPage
    from pairinteraction_gui.qobjects.item import ParameterItemBase


RangesKeys = Literal["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Distance", "Angle"]


class SystemConfig(BaseConfig):
    """Section for configuring the system parameters."""

    margin = (5, 15, 5, 5)
    spacing = 10

    title = "System"
    page: OneAtomPage | TwoAtomsPage | C6Page

    def setupEField(self, *, single_parameter: bool = False) -> None:
        efield_label = QLabel("<b>Electric field</b>")
        self.layout().addWidget(efield_label)

        parameter_class = ParameterItem if single_parameter else ParameterItemRange
        self.Ex = parameter_class(self, "Ex", unit="V/cm", tooltip_label="electric field in x-direction")
        self.Ey = parameter_class(self, "Ey", unit="V/cm", tooltip_label="electric field in y-direction")
        self.Ez = parameter_class(self, "Ez", unit="V/cm", tooltip_label="electric field in z-direction")

        self.layout().addWidget(self.Ex)
        self.layout().addWidget(self.Ey)
        self.layout().addWidget(self.Ez)

    def setupBField(self, *, single_parameter: bool = False) -> None:
        bfield_label = QLabel("<b>Magnetic field</b>")
        self.layout().addWidget(bfield_label)

        parameter_class = ParameterItem if single_parameter else ParameterItemRange
        self.Bx = parameter_class(self, "Bx", unit="Gauss", tooltip_label="magnetic field in x-direction")
        self.By = parameter_class(self, "By", unit="Gauss", tooltip_label="magnetic field in y-direction")
        self.Bz = parameter_class(self, "Bz", unit="Gauss", tooltip_label="magnetic field in z-direction")

        self.layout().addWidget(self.Bx)
        self.layout().addWidget(self.By)
        self.layout().addWidget(self.Bz)

    def setupDiamagnetism(self) -> None:
        self.layout().addWidget(QLabel("<b>Diamagnetism</b>"))
        self.diamagnetism = CheckboxItem(self, "Enable diamagnetism", checked=True)
        self.layout().addWidget(self.diamagnetism)

    def get_parameter_dict(self) -> dict[RangesKeys, list[float]]:
        """Return the electric and magnetic field parameters."""
        steps = self.page.calculation_config.steps.value()
        all_parameters = self._get_all_parameters()
        ranges_min_max: dict[str, tuple[float, float]] = {}
        for item in all_parameters:
            if not item.isChecked():
                continue
            if isinstance(item, ParameterItemRange):
                ranges_min_max[item.label.text()] = item.values()
            elif isinstance(item, ParameterItem):
                ranges_min_max[item.label.text()] = (item.value(0), item.value(0))

        if len(ranges_min_max) == 0:
            ranges_min_max["Bz"] = (0, 0)

        return {key: np.linspace(value[0], value[1], steps).tolist() for key, value in ranges_min_max.items()}  # type: ignore [misc]

    def _get_all_parameters(self) -> list[ParameterItemBase]:
        """Return all range items."""
        return [self.Ex, self.Ey, self.Ez, self.Bx, self.By, self.Bz]


class SystemConfigOneAtom(SystemConfig):
    page: OneAtomPage

    def setupWidget(self) -> None:
        self.setupEField()
        self.setupBField()
        self.setupDiamagnetism()


class SystemConfigTwoAtoms(SystemConfig):
    page: TwoAtomsPage | C6Page

    def setupWidget(self) -> None:
        self.setupEField()
        self.setupBField()
        self.setupDiamagnetism()
        self.setupDistance()
        self.setupAngle()
        self.setupOrder()

    def setupDistance(self, *, single_parameter: bool = False) -> None:
        label = QLabel("<b>Distance</b>")
        self.layout().addWidget(label)

        self.distance: ParameterItem | ParameterItemRange
        if single_parameter:
            self.distance = ParameterItem(
                self, "Distance", vdefault=10, vrange=(0, np.inf), unit="<span>&mu;m</span>", checked=False
            )
        else:
            self.distance = ParameterItemRange(
                self, "Distance", vdefaults=(3, 8), vrange=(0, np.inf), unit="<span>&mu;m</span>"
            )
        self.layout().addWidget(self.distance)

    def setupAngle(self, *, single_parameter: bool = False) -> None:
        label = QLabel("<b>Angle</b> (0° = z-axis, 90° = x-axis)")
        self.layout().addWidget(label)

        self.angle: ParameterItem | ParameterItemRange
        if single_parameter:
            self.angle = ParameterItem(self, "Angle", vdefault=0, vrange=(0, 360), unit="degree", checkable=False)
        else:
            self.angle = ParameterItemRange(
                self, "Angle", vdefaults=(0, 0), vrange=(0, 360), unit="degree", checkable=False
            )

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

    def _get_all_parameters(self) -> list[ParameterItemBase]:
        """Return all range items."""
        return [*super()._get_all_parameters(), self.distance, self.angle]


class SystemConfigC6(SystemConfigTwoAtoms):
    page: C6Page

    def setupWidget(self) -> None:
        self.setupEField(single_parameter=True)
        self.setupBField(single_parameter=True)
        self.setupDiamagnetism()
        self.setupAngle(single_parameter=True)
        self.setupOrder()

    def _get_all_parameters(self) -> list[ParameterItemBase]:
        """Return all range items."""
        return [*SystemConfig._get_all_parameters(self), self.angle]
