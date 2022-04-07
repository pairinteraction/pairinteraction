# Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
#
# This file is part of the pairinteraction GUI.
#
# The pairinteraction GUI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The pairinteraction GUI is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the pairinteraction GUI. If not, see <http://www.gnu.org/licenses/>.

from PyQt5 import QtCore, QtGui
import locale
from abc import ABCMeta, abstractmethod
from .unitmanagement import Quantity
import collections.abc


# === Dictionary to manage the elements of the GUI ===

class GuiDict(collections.abc.MutableMapping, metaclass=ABCMeta):

    def __init__(self, ui):
        self.store = dict()
        self._setup(self.store, ui)

    @abstractmethod
    def _setup(self, store, ui):
        pass

    def __getattr__(self, key):
        return self.__getitem__(key)

    def __getitem__(self, key):
        widget = self.store[key]['widget']
        unit = self.store[key]['unit']

        value = None
        if isinstance(widget, QtGui.QComboBox):
            value = str(widget.currentText())
        elif isinstance(widget, QtGui.QSpinBox):
            value = int(widget.value())
        elif isinstance(widget, QtGui.QDoubleSpinBox):
            value = float(widget.value())
        elif isinstance(widget, QtGui.QLineEdit):
            try:
                value = locale.atof(str(widget.text()))
            except ValueError:
                value = None
        elif isinstance(widget, QtGui.QCheckBox):
            value = widget.checkState() == QtCore.Qt.Checked
        elif isinstance(widget, QtGui.QRadioButton):
            value = widget.isChecked()
        elif isinstance(widget, QtGui.QGroupBox):
            value = widget.isChecked()

        return Quantity(value, unit)

    def __setitem__(self, key, value):
        if not isinstance(value, Quantity):
            raise Exception("value has to be of type quantity")

        widget = self.store[key]['widget']
        unit = self.store[key]['unit']

        value = value.toUU().magnitude

        if isinstance(widget, QtGui.QComboBox):
            index = widget.findText(value)
            if index >= 0:
                widget.setCurrentIndex(index)
        elif isinstance(widget, QtGui.QSpinBox):
            widget.setValue(value)
        elif isinstance(widget, QtGui.QDoubleSpinBox):
            widget.setValue(value)
        elif isinstance(widget, QtGui.QLineEdit):
            if value is None:
                widget.setText("None")
            else:
                widget.setText(locale.str(value))
        elif isinstance(widget, QtGui.QCheckBox):
            if value:
                widget.setCheckState(QtCore.Qt.Checked)
            else:
                widget.setCheckState(QtCore.Qt.Unchecked)
        elif isinstance(widget, QtGui.QRadioButton):
            widget.setChecked(value)
        elif isinstance(widget, QtGui.QGroupBox):
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
        if s == 'None':
            return (QtGui.QValidator.Acceptable, s, pos)

        lastpos = -1
        for c in s.lower():
            try:
                lastpos = 'none'[lastpos + 1:].index(c)
            except ValueError:
                return super().validate(s, pos)

        return (QtGui.QValidator.Intermediate, s, pos)

    def fixup(self, s):
        return 'None'


class DoublepositiveValidator(QtGui.QDoubleValidator):

    def __init__(self, parent=None):
        super().__init__(parent)

    def validate(self, s, pos):
        status = super().validate(s, pos)

        if status[0] == QtGui.QValidator.Intermediate and len(s) > 0 and s[0] == '-':
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
