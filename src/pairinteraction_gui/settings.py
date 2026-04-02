# SPDX-FileCopyrightText: 2026 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, ClassVar

from PySide6.QtCore import QObject, QSettings
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QRadioButton,
    QSpinBox,
    QStackedWidget,
    QWidget,
)

from pairinteraction import _backend
from pairinteraction_gui.config.base_config import BaseConfig

if TYPE_CHECKING:
    from pathlib import Path


logger = logging.getLogger(__name__)


class SettingsManager(QObject):
    """Settings manager."""

    widget_mappers: ClassVar[dict[type, tuple[str, str, type]]] = {
        QCheckBox: ("isChecked", "setChecked", bool),
        QSpinBox: ("value", "setValue", int),
        QDoubleSpinBox: ("value", "setValue", float),
        QRadioButton: ("isChecked", "setChecked", bool),
        QComboBox: ("currentText", "setCurrentText", str),
    }

    def __init__(self, cache_dir: Path | None = None) -> None:
        super().__init__()
        if cache_dir is None:
            cache_dir = _backend.get_cache_directory()
        path = cache_dir / "gui_settings.ini"
        path.parent.mkdir(parents=True, exist_ok=True)
        self.settings = QSettings(str(path), QSettings.Format.IniFormat)

    def value(self, key: str, default: object = None, value_type: type | None = None) -> object:
        if value_type is not None:
            return self.settings.value(key, defaultValue=default, type=value_type)
        return self.settings.value(key, defaultValue=default)

    def set_value(self, key: str, value: object) -> None:
        self.settings.setValue(key, value)

    def _get_mapper(self, widget: QWidget) -> tuple[str, str, type] | tuple[None, None, None]:
        return next((m for c, m in self.widget_mappers.items() if isinstance(widget, c)), (None, None, None))

    def update_widgets_from_settings(self, widget_map: dict[str, QWidget], *, combos_only: bool = False) -> None:
        """Set widget states from stored settings values."""
        for name, widget in widget_map.items():
            if combos_only and not isinstance(widget, QComboBox):
                continue

            getter, setter, dtype = self._get_mapper(widget)
            if not getter:
                continue

            value = getattr(widget, getter)()
            stored = self.settings.value(name, value, type=dtype)
            if stored is None:
                continue

            if setter:
                try:
                    getattr(widget, setter)(stored)
                except Exception as e:
                    logger.warning("Failed to restore setting '%s' with value '%s': %s", name, stored, e)

    def update_settings_from_widgets(self, widget_map: dict[str, QWidget]) -> None:
        """Save widget states into settings."""
        for name, widget in widget_map.items():
            getter, _setter, _dtype = self._get_mapper(widget)
            if getter:
                value = getattr(widget, getter)()
                if value is not None:
                    self.settings.setValue(name, value)

    def save_widget_state(self, root: QWidget, group: str) -> None:
        """Write the current state of all named input widgets under `group`."""
        if not isinstance(root, BaseConfig):
            return
        widget_map = self.collect_widgets(root)
        self.settings.beginGroup(group)
        self.update_settings_from_widgets(widget_map)
        self.settings.endGroup()

    def restore_widget_state(self, root: QWidget, group: str) -> None:
        """Restore widget state (two-pass: combos first, then others)."""
        if not isinstance(root, BaseConfig):
            return
        widget_map = self.collect_widgets(root)
        self.settings.beginGroup(group)
        self.update_widgets_from_settings(widget_map, combos_only=True)
        widget_map = self.collect_widgets(root)
        self.update_widgets_from_settings(widget_map, combos_only=False)
        self.settings.endGroup()

    def collect_widgets(self, root: QWidget, widget_map: dict[str, QWidget] | None = None) -> dict[str, QWidget]:
        if widget_map is None:
            widget_map = {}
        for child in root.children():
            if not isinstance(child, QWidget):
                continue
            name = child.objectName()
            if name and any(isinstance(child, c) for c in self.widget_mappers):
                if name in widget_map:
                    logger.warning("Duplicate widget name '%s' found. Only the last one will be saved/restored.", name)
                widget_map[name] = child

            # For stacked widgets, only recurse into the currently shown page
            if isinstance(child, QStackedWidget):
                current = child.currentWidget()
                if current is not None:
                    self.collect_widgets(current, widget_map)
            else:
                self.collect_widgets(child, widget_map)

        return widget_map
