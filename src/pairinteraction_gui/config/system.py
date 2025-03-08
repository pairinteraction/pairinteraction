from typing import TYPE_CHECKING, Union

import numpy as np
from PySide6.QtWidgets import (
    QLabel,
)

from pairinteraction import (
    complex as pi_complex,
    real as pi_real,
)
from pairinteraction_gui.config.base_config import BaseConfig
from pairinteraction_gui.qobjects import DoubleSpinBox, IntSpinBox, RangeItem, WidgetForm, WidgetV, parse_html

if TYPE_CHECKING:
    from pairinteraction_gui.page import SystemAtomPage, SystemPairPage


class SystemBaseConfig(BaseConfig):
    """Section for configuring the system parameters."""

    margin = (5, 15, 5, 5)
    spacing = 10

    title = "System"
    page: Union["SystemAtomPage", "SystemPairPage"]

    def postSetupWidget(self) -> None:
        self.layout().addStretch()

    def setupEField(self) -> None:
        efield_label = QLabel("<b>Electric field</b>")
        efield_label.setWordWrap(True)
        self.layout().addWidget(efield_label)

        self.efield_range = RangEfield()
        self.layout().addWidget(self.efield_range)

    def setupBField(self) -> None:
        bfield_label = QLabel("<b>Magnetic field</b>")
        bfield_label.setWordWrap(True)
        self.layout().addWidget(bfield_label)

        self.bfield_range = RangeBfield()
        self.layout().addWidget(self.bfield_range)

    def setupSteps(self) -> None:
        # Steps configuration section
        steps_label = QLabel("<b>Calculation steps</b>")
        steps_label.setWordWrap(True)
        self.layout().addWidget(steps_label)

        steps_widget = WidgetForm()
        self.steps_spinbox = IntSpinBox(vmin=1, vmax=1000, vdefault=20, tooltip="Number of steps for the calculation")
        steps_widget.layout().addRow("Number of steps:", self.steps_spinbox)
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
            system = pi.SystemAtom(basis)  # type: ignore
            efield = [fields[key][step] for key in ["Ex", "Ey", "Ez"]]
            bfield = [fields[key][step] for key in ["Bx", "By", "Bz"]]
            system.set_electric_field(efield, unit="V/cm")
            system.set_magnetic_field(bfield, unit="G")
            systems.append(system)  # type: ignore

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
        self.setupSteps()

    def setupDistance(self) -> None:
        label = QLabel("<b>Distance</b>")
        label.setWordWrap(True)
        self.layout().addWidget(label)

        self.distance_range = RangeDistance()
        self.layout().addWidget(self.distance_range)

    def setupAngle(self) -> None:
        label = QLabel("<b>Angle</b> (0° = z-axis, 90° = x-axis)")
        label.setWordWrap(True)
        self.layout().addWidget(label)

        self.angle_range = RangeAngle()
        self.layout().addWidget(self.angle_range)

    def get_distance(self) -> "RangeObject":
        """Return the distance range."""
        steps = self.steps_spinbox.value()
        item = self.distance_range.items[0]
        if not item.isChecked():
            return RangeObject(np.inf, np.inf, steps)
        return RangeObject(*item.values(), steps)  # type: ignore

    def get_angle(self) -> "RangeObject":
        """Return the angle range."""
        steps = self.steps_spinbox.value()
        item = self.angle_range.items[0]
        if not item.isChecked():
            return RangeObject(0, 0, steps)
        return RangeObject(*item.values(), steps)  # type: ignore


class RangeBase(WidgetV):
    """Base class for range configuration."""

    margin = (10, 0, 10, 0)
    spacing = 5

    items: list[RangeItem]


class RangEfield(RangeBase):
    """Configuration for electric field ranges in x, y, and z directions."""

    default_deactivated = ["Ey"]

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

    default_deactivated = ["By"]

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
            self, vmin=0, vmax=999, vdefault=3, tooltip=f"Minimum distance in {parse_html('mu')}"
        )
        max_spinbox = DoubleSpinBox(
            self, vmin=0, vmax=999, vdefault=10, tooltip=f"Maximum distance in {parse_html('mu')}"
        )

        item = RangeItem(self, "R", min_spinbox, max_spinbox, unit=parse_html("mu"))
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
        return self.list[step]
