# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
from PySide6.QtWidgets import QLabel

from pairinteraction_gui.config.base_config import BaseConfig
from pairinteraction_gui.qobjects import Item, QnItemDouble, QnItemInt, RangeItem

if TYPE_CHECKING:
    from pairinteraction_gui.page import OneAtomPage, TwoAtomsPage


RangesKeys = Literal["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Distance", "Angle", "IonDistance", "IonAngle"]


class SystemConfig(BaseConfig):
    """Section for configuring the system parameters."""

    margin = (5, 15, 5, 5)
    spacing = 10

    title = "System"
    page: OneAtomPage | TwoAtomsPage

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

        self.diamagnetism = Item(self, "Enable diamagnetism", checked=True)
        self.layout().addWidget(self.diamagnetism)

    def setupIonInteraction(self) -> None:
        self.layout().addWidget(QLabel("<b>Ion interaction</b>"))

        self.ion_distance = RangeItem(
            self,
            "Distance",
            vdefaults=(3, 8),
            vrange=(0, np.inf),
            unit="<span>&mu;m</span>",
            checked=False,
            key="IonDistance",
        )
        self.layout().addWidget(self.ion_distance)

        self.ion_angle = RangeItem(
            self, "Angle", vdefaults=(0, 0), vrange=(0, 360), unit="degree", checked=False, key="IonAngle"
        )
        self.layout().addWidget(self.ion_angle)

        self.ion_order = QnItemInt(
            self,
            "Multipole order",
            vmin=2,
            vmax=3,
            vdefault=3,
            tooltip="Select the order of the ion interaction multipole expansion",
            checkable=False,
        )
        self.layout().addWidget(self.ion_order)

        self.ion_charge = QnItemDouble(self, "Charge", vmin=-10, vmax=10, vdefault=1, unit="e", checkable=False)
        self.layout().addWidget(self.ion_charge)

    def get_ranges_dict(self) -> dict[RangesKeys, list[float]]:
        """Return the electric and magnetic field ranges."""
        steps = self.page.calculation_config.steps.value()
        all_ranges = self._get_all_ranges()
        ranges_min_max: dict[str, tuple[float, float]] = {
            item.key: item.values() for item in all_ranges if item.isChecked()
        }
        if len(ranges_min_max) == 0:
            ranges_min_max["Bz"] = (0, 0)

        return {key: np.linspace(value[0], value[1], steps).tolist() for key, value in ranges_min_max.items()}  # type: ignore [misc]

    def _get_all_ranges(self) -> list[RangeItem]:
        """Return all range items."""
        return [self.Ex, self.Ey, self.Ez, self.Bx, self.By, self.Bz]


class SystemConfigOneAtom(SystemConfig):
    page: OneAtomPage

    def setupWidget(self) -> None:
        self.setupEField()
        self.setupBField()
        self.setupIonInteraction()

    def _get_all_ranges(self) -> list[RangeItem]:
        """Return all range items."""
        return [*super()._get_all_ranges(), self.ion_distance, self.ion_angle]


class SystemConfigTwoAtoms(SystemConfig):
    page: TwoAtomsPage

    def setupWidget(self) -> None:
        self.setupEField()
        self.setupBField()
        self.setupIonInteraction()
        self.setupRydbergInteraction()

    def setupRydbergInteraction(self) -> None:
        self.layout().addWidget(QLabel("<b>Rydberg interaction</b>"))

        self.distance = RangeItem(self, "Distance", vdefaults=(3, 8), vrange=(0, np.inf), unit="<span>&mu;m</span>")
        self.layout().addWidget(self.distance)

        self.angle = RangeItem(self, "Angle", vdefaults=(0, 0), vrange=(0, 360), unit="degree", checkable=False)
        self.layout().addWidget(self.angle)

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
        return [*super()._get_all_ranges(), self.ion_distance, self.ion_angle, self.distance, self.angle]
