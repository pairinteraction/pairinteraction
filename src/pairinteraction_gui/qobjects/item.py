# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import typing as t
from typing import TYPE_CHECKING, Callable, Generic, Optional, TypeVar, Union

import numpy as np
from PySide6.QtWidgets import (
    QCheckBox,
    QLabel,
    QSpacerItem,
    QWidget,
)

from pairinteraction_gui.qobjects.spin_boxes import DoubleSpinBox, HalfIntSpinBox, IntSpinBox
from pairinteraction_gui.qobjects.widget import WidgetH

if TYPE_CHECKING:
    P = TypeVar("P")

ValueType = TypeVar("ValueType", int, float, complex)


@t.runtime_checkable
class NotSet(t.Protocol):
    """Singleton for a not set value and type at the same time.

    See Also:
    https://stackoverflow.com/questions/77571796/how-to-create-singleton-object-which-could-be-used-both-as-type-and-value-simi

    """

    @staticmethod
    def __not_set() -> None: ...


class Item(WidgetH):
    margin = (20, 0, 20, 0)
    spacing = 10

    def __init__(
        self,
        parent: QWidget,
        label: str,
        checked: bool = True,
    ) -> None:
        self.checkbox = QCheckBox()
        self.checkbox.setChecked(checked)

        self.label = label
        self._label = QLabel(label)
        self._label.setMinimumWidth(25)

        super().__init__(parent)

    def setupWidget(self) -> None:
        self.layout().addWidget(self.checkbox)
        self.layout().addWidget(self._label)

    def postSetupWidget(self) -> None:
        self.layout().addStretch(1)

    def isChecked(self) -> bool:
        """Return the state of the checkbox."""
        return self.checkbox.isChecked()

    def setChecked(self, checked: bool) -> None:
        """Set the state of the checkbox."""
        self.checkbox.setChecked(checked)

    def connectAll(self, func: Callable[[], None]) -> None:
        """Connect the function to the spinbox.valueChanged signal."""
        self.checkbox.stateChanged.connect(lambda state: func())


class _QnItem(WidgetH, Generic[ValueType]):
    """Widget for displaying a range with min and max spinboxes."""

    margin = (20, 0, 20, 0)
    spacing = 10
    _spinbox_class: type[Union[IntSpinBox, HalfIntSpinBox, DoubleSpinBox]]

    def __init__(
        self,
        parent: QWidget,
        label: str,
        vmin: ValueType = 0,
        vmax: ValueType = 999,
        vdefault: ValueType = 0,
        vstep: Optional[ValueType] = None,
        unit: str = "",
        tooltip: Optional[str] = None,
        checkable: bool = True,
        checked: bool = True,
    ) -> None:
        tooltip = tooltip if tooltip is not None else f"{label} in {unit}"

        self.checkbox: Union[QCheckBox, QSpacerItem]
        if checkable:
            self.checkbox = QCheckBox()
            self.checkbox.setChecked(checked)
            self.checkbox.stateChanged.connect(self._on_checkbox_changed)
        else:
            self.checkbox = QSpacerItem(25, 0)

        self.label = QLabel(label)
        self.label.setMinimumWidth(25)

        self.spinbox = self._spinbox_class(parent, vmin, vmax, vdefault, vstep, tooltip=tooltip)  # type: ignore [arg-type]
        self.spinbox.setObjectName(f"{label.lower()}")
        self.spinbox.setMinimumWidth(100)

        self.unit = QLabel(unit)

        super().__init__(parent)

    def setupWidget(self) -> None:
        if isinstance(self.checkbox, QCheckBox):
            self.layout().addWidget(self.checkbox)
        elif isinstance(self.checkbox, QSpacerItem):
            self.layout().addItem(self.checkbox)
        self.layout().addWidget(self.label)
        self.layout().addWidget(self.spinbox)
        self.layout().addWidget(self.unit)

    def postSetupWidget(self) -> None:
        self.layout().addStretch(1)
        self._on_checkbox_changed(self.isChecked())

    def connectAll(self, func: Callable[[], None]) -> None:
        """Connect the function to the spinbox.valueChanged signal."""
        if isinstance(self.checkbox, QCheckBox):
            self.checkbox.stateChanged.connect(lambda state: func())
        self.spinbox.valueChanged.connect(lambda value: func())

    def _on_checkbox_changed(self, state: bool) -> None:
        """Update the enabled state of widget when checkbox changes."""
        self.spinbox.setEnabled(state)

    def isChecked(self) -> bool:
        """Return the state of the checkbox."""
        if isinstance(self.checkbox, QSpacerItem):
            return True
        return self.checkbox.isChecked()

    def setChecked(self, checked: bool) -> None:
        """Set the state of the checkbox."""
        if isinstance(self.checkbox, QSpacerItem):
            if checked:
                return
            raise ValueError("Cannot uncheck a non-checkable item.")
        self.checkbox.setChecked(checked)

    def value(
        self,
        default: Union["P", NotSet] = NotSet,
    ) -> Union[ValueType, "P"]:
        """Return the value of the spinbox."""
        if not self.isChecked():
            if isinstance(default, NotSet):
                raise ValueError("Checkbox is not checked and no default value is provided.")
            return default
        return self.spinbox.value()  # type: ignore [return-value]

    def setValue(self, value: ValueType) -> None:
        """Set the value of the spinbox and set the checkbox state to checked if applicable."""
        if isinstance(self.checkbox, QCheckBox):
            self.checkbox.setChecked(True)
        self.spinbox.setValue(value)  # type: ignore [arg-type]


class QnItemInt(_QnItem[int]):
    _spinbox_class = IntSpinBox


class QnItemHalfInt(_QnItem[float]):
    _spinbox_class = HalfIntSpinBox


class QnItemDouble(_QnItem[float]):
    _spinbox_class = DoubleSpinBox


class RangeItem(WidgetH):
    """Widget for displaying a range with min and max spinboxes."""

    margin = (20, 0, 20, 0)
    spacing = 10

    def __init__(
        self,
        parent: QWidget,
        label: str,
        vdefaults: tuple[float, float] = (0, 0),
        vrange: tuple[float, float] = (-np.inf, np.inf),
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

    def setChecked(self, checked: bool) -> None:
        """Set the state of the checkbox."""
        if isinstance(self.checkbox, QSpacerItem):
            if checked:
                return
            raise ValueError("Cannot uncheck a non-checkable item.")
        self.checkbox.setChecked(checked)

    def values(self, default: Union[tuple[float, float], "P", NotSet] = NotSet) -> Union[tuple[float, float], "P"]:
        """Return the values of the min and max spinboxes."""
        if not self.isChecked():
            if isinstance(default, NotSet):
                raise ValueError("Checkbox is not checked and no default value is provided.")
            return default
        return (self.min_spinbox.value(), self.max_spinbox.value())

    def setValues(self, vmin: float, vmax: float) -> None:
        """Set the values of the min and max spinboxes and set the checkbox state to checked if applicable."""
        if isinstance(self.checkbox, QCheckBox):
            self.checkbox.setChecked(True)
        self.min_spinbox.setValue(vmin)
        self.max_spinbox.setValue(vmax)
