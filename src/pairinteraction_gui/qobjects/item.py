from typing import Callable, Union

from PySide6.QtWidgets import (
    QCheckBox,
    QLabel,
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

        self.layout().addStretch()

        # Initial state
        self._on_checkbox_changed(self.isChecked())

    def _on_checkbox_changed(self, state: bool) -> None:
        """Update the enabled state of widgets when checkbox changes."""
        self.setStyleSheet("color: black" if state else "color: gray")
        for spinbox in self.spinboxes.values():
            spinbox.setEnabled(state)

    def isChecked(self) -> bool:
        """Return the state of the checkbox."""
        return self.checkbox.isChecked()

    def connectAll(self, func: Callable) -> None:
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


class RangeItem(Item):
    """Widget for displaying a range with min and max spinboxes."""

    def __init__(
        self,
        parent: QWidget,
        label: str,
        min_spinbox: Union[IntSpinBox, HalfIntSpinBox, DoubleSpinBox],
        max_spinbox: Union[IntSpinBox, HalfIntSpinBox, DoubleSpinBox],
        unit: str = "",
        checked: bool = True,
    ) -> None:
        min_spinbox.setObjectName(f"{label}_min")
        max_spinbox.setObjectName(f"{label}_max")
        spinboxes = {"min": min_spinbox, "max": max_spinbox}
        super().__init__(parent, label, spinboxes, unit, checked)

    def setupWidget(self) -> None:
        super().setupWidget()
        self.layout().addWidget(self.spinboxes["min"])
        _to_label = QLabel("to")
        self.layout().addWidget(_to_label)
        self.layout().addWidget(self.spinboxes["max"])

    def values(self) -> tuple[float, float]:
        """Return the values of the min and max spinboxes."""
        return self.spinboxes["min"].value(), self.spinboxes["max"].value()
