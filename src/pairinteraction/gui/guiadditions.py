# SPDX-FileCopyrightText: 2024 Sebastian Weber, Henri Menke, and contributors. All rights reserved.
# SPDX-License-Identifier: GPL-3.0-or-later
import collections.abc
import locale
from abc import ABCMeta, abstractmethod

from PyQt5 import QtCore, QtGui, QtWidgets

from pairinteraction.gui.unitmanagement import Quantity

# === Dictionary to manage the elements of the GUI ===


class GuiDict(collections.abc.MutableMapping, metaclass=ABCMeta):
    def __init__(self, ui):
        self.store = {}
        self._setup(self.store, ui)

    @abstractmethod
    def _setup(self, store, ui):
        pass

    def __getattr__(self, key):
        return self.__getitem__(key)

    def __getitem__(self, key):
        widget = self.store[key]["widget"]
        unit = self.store[key]["unit"]

        value = None
        if isinstance(widget, QtWidgets.QComboBox):
            value = str(widget.currentText())
        elif isinstance(widget, QtWidgets.QSpinBox):
            value = int(widget.value())
        elif isinstance(widget, QtWidgets.QDoubleSpinBox):
            value = float(widget.value())
        elif isinstance(widget, QtWidgets.QLineEdit):
            try:
                value = locale.atof(str(widget.text()))
            except ValueError:
                value = None
        elif isinstance(widget, QtWidgets.QCheckBox):
            value = widget.checkState() == QtCore.Qt.Checked
        elif isinstance(widget, QtWidgets.QRadioButton):
            value = widget.isChecked()
        elif isinstance(widget, QtWidgets.QGroupBox):
            value = widget.isChecked()

        return Quantity(value, unit)

    def __setitem__(self, key, value):
        if not isinstance(value, Quantity):
            raise Exception("value has to be of type quantity")

        widget = self.store[key]["widget"]

        value = value.toUU().magnitude

        if isinstance(widget, QtWidgets.QComboBox):
            index = widget.findText(value)
            if index >= 0:
                widget.setCurrentIndex(index)
        elif isinstance(widget, QtWidgets.QSpinBox):
            widget.setValue(value)
        elif isinstance(widget, QtWidgets.QDoubleSpinBox):
            widget.setValue(value)
        elif isinstance(widget, QtWidgets.QLineEdit):
            if value is None:
                widget.setText("None")
            else:
                widget.setText(locale.str(value))
        elif isinstance(widget, QtWidgets.QCheckBox):
            if value:
                widget.setCheckState(QtCore.Qt.Checked)
            else:
                widget.setCheckState(QtCore.Qt.Unchecked)
        elif isinstance(widget, QtWidgets.QRadioButton):
            widget.setChecked(value)
        elif isinstance(widget, QtWidgets.QGroupBox):
            widget.setChecked(value)

    def __delitem__(self, key):
        del self.store[key]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)


# === Validators ===


class DoublenoneValidator(QtGui.QDoubleValidator):
    def __init__(self, parent=None):
        super().__init__(parent)

    def validate(self, s, pos):
        if s == "None":
            return (QtGui.QValidator.Acceptable, s, pos)

        lastpos = -1
        for c in s.lower():
            try:
                lastpos = "none"[lastpos + 1 :].index(c)
            except ValueError:
                return super().validate(s, pos)

        return (QtGui.QValidator.Intermediate, s, pos)

    def fixup(self, s):
        return "None"


class DoublepositiveValidator(QtGui.QDoubleValidator):
    def __init__(self, parent=None):
        super().__init__(parent)

    def validate(self, s, pos):
        status = super().validate(s, pos)

        if status[0] == QtGui.QValidator.Intermediate and len(s) > 0 and s[0] == "-":
            return (QtGui.QValidator.Invalid, s, pos)

        if status[0] == QtGui.QValidator.Acceptable and locale.atof(s) < 0:
            return (QtGui.QValidator.Invalid, s, pos)

        return status

    def fixup(self, s):
        return "0"


class DoubledeltaValidator(QtGui.QDoubleValidator):
    def __init__(self, parent=None):
        super().__init__(parent)

    def validate(self, s, pos):
        status = super().validate(s, pos)

        if status[0] == QtGui.QValidator.Acceptable and locale.atof(s) < 0 and locale.atof(s) != -1:
            return (QtGui.QValidator.Intermediate, s, pos)

        return status

    def fixup(self, s):
        if locale.atof(s) < 0:
            return "-1"
        return "0"


class DoubleValidator(QtGui.QDoubleValidator):
    def __init__(self, parent=None):
        super().__init__(parent)

    def validate(self, s, pos):
        return super().validate(s, pos)

    def fixup(self, s):
        return "0"
