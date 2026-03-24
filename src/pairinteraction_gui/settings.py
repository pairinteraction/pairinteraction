# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import json
import logging
from typing import TYPE_CHECKING, Any

from PySide6.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox, QSpinBox, QStackedWidget, QWidget

from pairinteraction import _backend
from pairinteraction_gui.config.base_config import BaseConfig

if TYPE_CHECKING:
    from pathlib import Path


logger = logging.getLogger(__name__)


def get_settings_path() -> Path:
    """Return the path to the GUI settings file."""
    return _backend.get_cache_directory() / "gui_settings.json"


def load_settings() -> dict[str, Any]:
    """Load GUI settings from disk, returning an empty dict on any error."""
    try:
        with get_settings_path().open() as f:
            return json.load(f)  # type: ignore[no-any-return]
    except (FileNotFoundError, json.JSONDecodeError, OSError) as e:
        logger.warning("Could not load GUI settings: %s", e)
        return {}


def save_settings(data: dict[str, Any]) -> None:
    """Save GUI settings to disk."""
    path = get_settings_path()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as f:
            json.dump(data, f, indent=2)
    except OSError as e:
        logger.warning("Could not save GUI settings: %s", e)


def collect_widget_state(root: QWidget) -> dict[str, Any]:
    """Collect {objectName: value} for all named input widgets.

    Only enters the current page of QStackedWidgets to avoid collecting
    hidden pages (e.g. non-selected species widgets).
    """
    if not isinstance(root, BaseConfig):
        return {}
    state: dict[str, Any] = {}
    _collect(root, state)
    return state


def _collect(widget: QWidget, state: dict[str, Any]) -> None:
    for child in widget.children():
        if not isinstance(child, QWidget):
            continue
        name = child.objectName()
        if name:
            if isinstance(child, QComboBox):
                state[name] = child.currentText()
            elif isinstance(child, (QSpinBox, QDoubleSpinBox)):
                state[name] = child.value()
            elif isinstance(child, QCheckBox):
                state[name] = child.isChecked()
        # For stacked widgets, only recurse into the currently shown page
        if isinstance(child, QStackedWidget):
            current = child.currentWidget()
            if current is not None:
                _collect(current, state)
        else:
            _collect(child, state)


def restore_widget_state(root: QWidget, state: dict[str, Any]) -> None:
    """Restore widget state from a {objectName: value} dict.

    Two-pass: first restore QComboBoxes (which may create new stacked pages
    via species-change signals), then restore spinboxes and checkboxes.
    """
    if len(state) == 0 or not isinstance(root, BaseConfig):
        return
    _restore(root, state, combos_only=True)
    _restore(root, state, combos_only=False)


def _restore(widget: QWidget, state: dict[str, Any], *, combos_only: bool) -> None:  # noqa: C901
    for child in widget.children():
        if not isinstance(child, QWidget):
            continue
        name = child.objectName()
        if name and name in state:
            if isinstance(child, QComboBox):
                index = child.findText(state[name])
                if index != -1:
                    child.setCurrentIndex(index)
                else:
                    logger.warning("Could not restore setting '%s': '%s' not found in combo box", name, state[name])
            elif not combos_only:
                if isinstance(child, (QSpinBox, QDoubleSpinBox)):
                    child.setValue(state[name])
                elif isinstance(child, QCheckBox):
                    child.setChecked(state[name])
        if isinstance(child, QStackedWidget):
            current = child.currentWidget()
            if current is not None:
                _restore(current, state, combos_only=combos_only)
        else:
            _restore(child, state, combos_only=combos_only)
