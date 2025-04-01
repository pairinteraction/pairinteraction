# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, ClassVar, Union

import numpy as np
from PySide6.QtWidgets import (
    QLabel,
)

from pairinteraction import (
    complex as pi_complex,
    real as pi_real,
)
from pairinteraction_gui.config.base_config import BaseConfig
from pairinteraction_gui.qobjects import DoubleSpinBox, IntSpinBox, RangeItem, WidgetForm, WidgetV

if TYPE_CHECKING:
    from pairinteraction_gui.page import OneAtomPage, TwoAtomsPage


class SystemBaseConfig(BaseConfig):
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

        self.efield_range = RangEfield()
        self.layout().addWidget(self.efield_range)

    def setupBField(self) -> None:
        bfield_label = QLabel("<b>Magnetic field</b>")
        self.layout().addWidget(bfield_label)

        self.bfield_range = RangeBfield()
        self.layout().addWidget(self.bfield_range)

    def setupSteps(self) -> None:
        steps_label = QLabel("<b>Calculation steps</b>")
        self.layout().addWidget(steps_label)

        steps_widget = WidgetForm(margin=(0, 0, 0, 0))
        self.steps_spinbox = IntSpinBox(vmin=1, vmax=10000, vdefault=20, tooltip="Number of steps for the calculation")
        steps_widget.layout().addRow("Number of steps", self.steps_spinbox)
        self.layout().addWidget(steps_widget)

    def get_fields(self) -> dict[str, "RangeObject"]:
        """Return the electric and magnetic field ranges."""
        steps = self.steps_spinbox.value()
        all_fields: dict[str, tuple[float, float]] = {
            item.label: item.values() if item.isChecked() else (0, 0)
            for item in self.efield_range.items + self.bfield_range.items
        }
        fields: dict[str, RangeObject] = {}
        for key, value in all_fields.items():
            fields[key] = RangeObject(*value, steps)
        return fields

    def get_systems(self, atom: int) -> Union[list["pi_real.SystemAtom"], list["pi_complex.SystemAtom"]]:
        """Return a list of system atoms and corresponding field min/max values for the given configuration."""
        fields = self.get_fields()
        steps = next(iter(fields.values())).steps

        isreal = all(f.is_zero() for key, f in fields.items() if key.endswith("y"))
        pi = pi_real if isreal else pi_complex

        basis = self.page.basis_config.get_basis(atom, "real" if isreal else "complex")
        systems: list[pi_real.SystemAtom] = []
        for step in range(steps):
            system = pi.SystemAtom(basis)
            efield = [fields[key][step] for key in ["Ex", "Ey", "Ez"]]
            bfield = [fields[key][step] for key in ["Bx", "By", "Bz"]]
            system.set_electric_field(efield, unit="V/cm")
            system.set_magnetic_field(bfield, unit="G")
            systems.append(system)

        return systems


class SystemAtomConfig(SystemBaseConfig):
    def setupWidget(self) -> None:
        self.setupEField()
        self.setupBField()
        self.setupSteps()


class SystemPairConfig(SystemBaseConfig):
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

        self.distance_range = RangeDistance()
        self.layout().addWidget(self.distance_range)

    def setupAngle(self) -> None:
        label = QLabel("<b>Angle</b> (0° = z-axis, 90° = x-axis)")
        self.layout().addWidget(label)

        self.angle_range = RangeAngle()
        self.layout().addWidget(self.angle_range)

    def setupOrder(self) -> None:
        steps_label = QLabel("<b>Multipole expansion order</b>")
        self.layout().addWidget(steps_label)

        steps_widget = WidgetForm(margin=(0, 0, 0, 0))
        self.order = IntSpinBox(self, vmin=3, vmax=5, vdefault=3, tooltip="Select the order of the multipole expansion")
        steps_widget.layout().addRow("Multipole order", self.order)
        self.layout().addWidget(steps_widget)

    def get_distance(self) -> "RangeObject":
        """Return the distance range."""
        steps = self.steps_spinbox.value()
        item = self.distance_range.items[0]
        if not item.isChecked():
            return RangeObject(np.inf, np.inf, steps)
        values = item.values()
        return RangeObject(values[0], values[1], steps)

    def get_angle(self) -> "RangeObject":
        """Return the angle range."""
        steps = self.steps_spinbox.value()
        item = self.angle_range.items[0]
        if not item.isChecked():
            return RangeObject(0, 0, steps)
        values = item.values()
        return RangeObject(values[0], values[1], steps)

    def get_order(self) -> int:
        """Return the order of the multipole expansion."""
        return int(self.order.value())


class RangeBase(WidgetV):
    """Base class for range configuration."""

    margin = (10, 0, 10, 0)
    spacing = 5

    items: list[RangeItem]
    default_deactivated: ClassVar[list[str]] = []


class RangEfield(RangeBase):
    """Configuration for electric field ranges in x, y, and z directions."""

    default_deactivated: ClassVar = ["Ey"]

    def setupWidget(self) -> None:
        self.items = []

        for x in "xyz":
            min_spinbox = DoubleSpinBox(
                self, vmin=-999, vmax=999, tooltip=f"Minimum electric field in {x}-direction (V/cm)"
            )
            max_spinbox = DoubleSpinBox(
                self, vmin=-999, vmax=999, tooltip=f"Maximum electric field in {x}-direction (V/cm)"
            )

            checked = f"E{x}" not in self.default_deactivated
            item = RangeItem(self, f"E{x}", min_spinbox, max_spinbox, "V/cm", checked)
            self.items.append(item)
            self.layout().addWidget(item)


class RangeBfield(RangeBase):
    """Configuration for magnetic field ranges in x, y, and z directions."""

    default_deactivated: ClassVar = ["By"]

    def setupWidget(self) -> None:
        self.items = []

        for x in "xyz":
            min_spinbox = DoubleSpinBox(
                self, vmin=-999, vmax=999, tooltip=f"Minimum magnetic field in {x}-direction (Gauss)"
            )
            max_spinbox = DoubleSpinBox(
                self, vmin=-999, vmax=999, tooltip=f"Maximum magnetic field in {x}-direction (Gauss)"
            )

            checked = f"B{x}" not in self.default_deactivated
            item = RangeItem(self, f"B{x}", min_spinbox, max_spinbox, "Gauss", checked)
            self.items.append(item)
            self.layout().addWidget(item)


class RangeDistance(RangeBase):
    """Configuration for the distance."""

    def setupWidget(self) -> None:
        self.items = []

        min_spinbox = DoubleSpinBox(
            self, vmin=0, vmax=999, vdefault=3, tooltip="Minimum distance in <span>&mu;m</span>"
        )
        max_spinbox = DoubleSpinBox(
            self, vmin=0, vmax=999, vdefault=10, tooltip="Maximum distance in <span>&mu;m</span>"
        )

        item = RangeItem(self, "R", min_spinbox, max_spinbox, unit="<span>&mu;m</span>")
        self.items.append(item)
        self.layout().addWidget(item)


class RangeAngle(RangeBase):
    """Configuration for the angle."""

    def setupWidget(self) -> None:
        self.items = []

        min_spinbox = DoubleSpinBox(self, vmin=0, vmax=360, vdefault=0, tooltip="Minimum angle in °")
        max_spinbox = DoubleSpinBox(self, vmin=0, vmax=360, vdefault=0, tooltip="Maximum angle in °")

        item = RangeItem(self, "Angle", min_spinbox, max_spinbox, unit="degrees")
        self.items.append(item)
        self.layout().addWidget(item)


class RangeObject:
    def __init__(self, min_value: float, max_value: float, steps: int) -> None:
        self.min_value = min_value
        self.max_value = max_value
        self.steps = steps
        self.list = np.linspace(self.min_value, self.max_value, self.steps)

    def is_zero(self) -> bool:
        return self.min_value == 0 and self.max_value == 0

    def is_constant(self) -> bool:
        return self.min_value == self.max_value

    def __getitem__(self, step: int) -> float:
        return self.list[step]  # type: ignore [no-any-return] # FIXME, no idea why numpy type hints dont work here
