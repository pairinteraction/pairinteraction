from typing import Callable, Optional, Union

from PySide6.QtWidgets import QCheckBox, QLabel, QSpacerItem, QWidget

from pairinteraction_gui.qobjects.spin_boxes import DoubleSpinBox
from pairinteraction_gui.qobjects.widget import WidgetH


class RangeItem(WidgetH):
    """Widget for displaying a range with min and max spinboxes."""

    margin = (20, 0, 20, 0)
    spacing = 10

    def __init__(
        self,
        parent: QWidget,
        short_label: str,
        value_range: tuple[float, float] = (0, 999),
        value_defaults: tuple[float, float] = (0, 0),
        unit: str = "",
        long_label: Optional[str] = None,
        checkable: bool = True,
        checked: bool = True,
    ) -> None:
        long_label = long_label if long_label is not None else short_label

        self.checkbox: Union[QCheckBox, QSpacerItem]
        if checkable:
            self.checkbox = QCheckBox()
            self.checkbox.setChecked(checked)
            self.checkbox.stateChanged.connect(self._on_checkbox_changed)
        else:
            self.checkbox = QSpacerItem(0, 0)

        self.label = QLabel(short_label)
        self.label.setMinimumWidth(25)

        self.min_spinbox = DoubleSpinBox(
            parent, *value_range, value_defaults[0], tooltip=f"Minimum {long_label} in {unit}"
        )
        self.max_spinbox = DoubleSpinBox(
            parent, *value_range, value_defaults[1], tooltip=f"Maximum {long_label} in {unit}"
        )
        self.min_spinbox.setObjectName(f"{short_label.lower()}_min")
        self.max_spinbox.setObjectName(f"{short_label.lower()}_max")

        self.unit = QLabel(unit)

        super().__init__(parent)

    def setupWidget(self) -> None:
        self.layout().addWidget(self.checkbox)
        self.layout().addWidget(self.label)

        self.layout().addWidget(self.min_spinbox)
        self.layout().addWidget(QLabel("to"))
        self.layout().addWidget(self.max_spinbox)

        self.layout().addWidget(self.unit)
        self.layout().addStretch()

    def postSetupWidget(self) -> None:
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
        self.setStyleSheet("color: black" if state else "color: gray")
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
