# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import Callable, Optional, Union

from PySide6.QtWidgets import (
    QCheckBox,
    QLabel,
    QSpacerItem,
    QWidget,
)

from pairinteraction_gui.qobjects import parse_html
from pairinteraction_gui.qobjects.spin_boxes import DoubleSpinBox, HalfIntSpinBox, IntSpinBox
from pairinteraction_gui.qobjects.widget import WidgetH


class Item(WidgetH):
    margin = (20, 0, 20, 0)
    spacing = 10

    def __init__(
        self,
        parent: QWidget,
        label: str,
        spinboxes: dict[str, Union[IntSpinBox, HalfIntSpinBox, DoubleSpinBox]],
        unit: str = "",
        checked: bool = True,
    ) -> None:
        self.checkbox = QCheckBox()
        self.checkbox.setChecked(checked)
        self.checkbox.stateChanged.connect(self._on_checkbox_changed)

        self.label = label
        self._label = QLabel(parse_html(label))
        self._label.setMinimumWidth(25)

        self.spinboxes = spinboxes

        self.unit = unit

        super().__init__(parent)

    def setupWidget(self) -> None:
        self.layout().addWidget(self.checkbox)
        self.layout().addWidget(self._label)

    def postSetupWidget(self) -> None:
        self._unit = None
        if self.unit is not None and len(self.unit) > 0:
            self._unit = QLabel(self.unit)
            self.layout().addWidget(self._unit)

        self.layout().addStretch(1)

        # Initial state
        self._on_checkbox_changed(self.isChecked())

    def _on_checkbox_changed(self, state: bool) -> None:
        """Update the enabled state of widgets when checkbox changes."""
        for spinbox in self.spinboxes.values():
            spinbox.setEnabled(state)

    def isChecked(self) -> bool:
        """Return the state of the checkbox."""
        return self.checkbox.isChecked()

    def connectAll(self, func: Callable[[], None]) -> None:
        """Connect the function to the spinbox.valueChanged signal."""
        self.checkbox.stateChanged.connect(lambda state: func())
        for spinbox in self.spinboxes.values():
            spinbox.valueChanged.connect(lambda value: func())


class QnItem(Item):
    """Widget for displaying a Qn value with a single spinbox."""

    def __init__(
        self,
        parent: QWidget,
        label: str,
        spinbox: Union[IntSpinBox, HalfIntSpinBox, DoubleSpinBox],
        unit: str = "",
        checked: bool = True,
    ) -> None:
        spinbox.setObjectName(label)
        spinbox.setMinimumWidth(100)

        spinboxes = {"value": spinbox}
        super().__init__(parent, label, spinboxes, unit, checked)

    def setupWidget(self) -> None:
        super().setupWidget()
        self.layout().addWidget(self.spinboxes["value"])

    def value(self) -> Union[int, float]:
        """Return the value of the spinbox."""
        return self.spinboxes["value"].value()


class RangeItem(WidgetH):
    """Widget for displaying a range with min and max spinboxes."""

    margin = (20, 0, 20, 0)
    spacing = 10

    def __init__(
        self,
        parent: QWidget,
        label: str,
        vdefaults: tuple[float, float] = (0, 0),
        vrange: tuple[float, float] = (-999, 999),
        unit: str = "",
        tooltip_label: Optional[str] = None,
        checkable: bool = True,
        checked: bool = True,
    ) -> None:
        tooltip_label = tooltip_label if tooltip_label is not None else label

        self.checkbox: Union[QCheckBox, QSpacerItem]
        if checkable:
            self.checkbox = QCheckBox()
            self.checkbox.setChecked(checked)
            self.checkbox.stateChanged.connect(self._on_checkbox_changed)
        else:
            self.checkbox = QSpacerItem(25, 0)

        self.label = QLabel(label)
        self.label.setMinimumWidth(25)

        self.min_spinbox = DoubleSpinBox(parent, *vrange, vdefaults[0], tooltip=f"Minimum {tooltip_label} in {unit}")
        self.max_spinbox = DoubleSpinBox(parent, *vrange, vdefaults[1], tooltip=f"Maximum {tooltip_label} in {unit}")
        self.min_spinbox.setObjectName(f"{label.lower()}_min")
        self.max_spinbox.setObjectName(f"{label.lower()}_max")

        self.unit = QLabel(unit)

        super().__init__(parent)

    def setupWidget(self) -> None:
        if isinstance(self.checkbox, QCheckBox):
            self.layout().addWidget(self.checkbox)
        elif isinstance(self.checkbox, QSpacerItem):
            self.layout().addItem(self.checkbox)
        self.layout().addWidget(self.label)

        self.layout().addWidget(self.min_spinbox)
        self.layout().addWidget(QLabel("to"))
        self.layout().addWidget(self.max_spinbox)

        self.layout().addWidget(self.unit)

    def postSetupWidget(self) -> None:
        self.layout().addStretch(1)
        self._on_checkbox_changed(self.isChecked())

    @property
    def spinboxes(self) -> tuple[DoubleSpinBox, DoubleSpinBox]:
        """Return the min and max spinboxes."""
        return (self.min_spinbox, self.max_spinbox)

    def connectAll(self, func: Callable[[], None]) -> None:
        """Connect the function to the spinbox.valueChanged signal."""
        if isinstance(self.checkbox, QCheckBox):
            self.checkbox.stateChanged.connect(lambda state: func())
        for spinbox in self.spinboxes:
            spinbox.valueChanged.connect(lambda value: func())

    def _on_checkbox_changed(self, state: bool) -> None:
        """Update the enabled state of widgets when checkbox changes."""
        for spinbox in self.spinboxes:
            spinbox.setEnabled(state)

    def isChecked(self) -> bool:
        """Return the state of the checkbox."""
        if isinstance(self.checkbox, QSpacerItem):
            return True
        return self.checkbox.isChecked()

    def values(self) -> tuple[float, float]:
        """Return the values of the min and max spinboxes."""
        return (self.min_spinbox.value(), self.max_spinbox.value())
