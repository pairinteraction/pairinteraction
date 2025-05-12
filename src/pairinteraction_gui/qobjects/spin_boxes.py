# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import Optional

from PySide6.QtWidgets import (
    QDoubleSpinBox,
    QSpinBox,
    QWidget,
)


class IntSpinBox(QSpinBox):
    """Custom spin box for integer values."""

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        vmin: int = 0,
        vmax: int = 999,
        vdefault: int = 0,
        vstep: Optional[int] = None,
        suffix: Optional[str] = None,
        tooltip: Optional[str] = None,
    ) -> None:
        """Initialize the integer spin box."""
        super().__init__(parent)

        self.setRange(int(vmin), int(vmax))
        self.setValue(int(vdefault))
        self.setSingleStep(vstep if vstep is not None else 1)

        if suffix:
            self.setSuffix(suffix)
        if tooltip:
            self.setToolTip(tooltip)


class HalfIntSpinBox(QDoubleSpinBox):
    """Custom spin box for half int values."""

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        vmin: float = 0.5,
        vmax: float = 999.5,
        vdefault: float = 0.5,
        vstep: Optional[int] = None,
        suffix: Optional[str] = None,
        tooltip: Optional[str] = None,
    ) -> None:
        """Initialize the double spin box."""
        super().__init__(parent)

        vstep = vstep if vstep is not None else 1
        assert vdefault % 1 == 0.5, "Default value must be a half integer."  # NOSONAR
        assert vstep % 1 == 0, "Step value must be an integer."

        self.setRange(vmin, vmax)
        self.setValue(vdefault)
        self.setSingleStep(vstep)
        self.setDecimals(1)

        if suffix:
            self.setSuffix(suffix)
        if tooltip:
            self.setToolTip(tooltip)

    def valueFromText(self, text: str) -> float:
        """Convert text to value, ensuring it's a half integer."""
        value = super().valueFromText(text)
        return round(value - 0.49) + 0.5


class DoubleSpinBox(QDoubleSpinBox):
    """Custom spin box for double values."""

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        vmin: float = 0,
        vmax: float = 999.9,
        vdefault: float = 0,
        vstep: Optional[float] = None,
        suffix: Optional[str] = None,
        decimals: int = 1,
        tooltip: Optional[str] = None,
    ) -> None:
        """Initialize the double spin box."""
        super().__init__(parent)

        self.setDecimals(decimals)
        self.setSingleStep(vstep if vstep is not None else 10**-decimals)
        self.setRange(vmin, vmax)
        self.setValue(vdefault)

        if suffix:
            self.setSuffix(suffix)
        if tooltip:
            self.setToolTip(tooltip)
