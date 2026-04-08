# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import importlib.util
import logging
import sys
from pathlib import Path
from typing import TYPE_CHECKING

from PySide6.QtCore import QFileSystemWatcher, QObject, QTimer, Signal

from pairinteraction_gui.theme.palette import build_application_palette

if TYPE_CHECKING:
    from types import ModuleType

    from PySide6.QtGui import QPalette

logger = logging.getLogger(__name__)

_THEME_FILE_NAME = "theme.qss"
_PALETTE_FILE_NAME = "palette.py"


class ThemeSignals(QObject):
    """Signals emitted when theme files change."""

    themes_reloaded = Signal()


class ThemeManager(QObject):
    """Load theme files and optionally watch them for development reloads."""

    def __init__(self, theme_dir: Path | None = None) -> None:
        super().__init__()
        self.theme_dir = theme_dir or Path(__file__).parent
        self._theme = (self.theme_dir / _THEME_FILE_NAME).read_text()
        self._palette = build_application_palette()

        # Variables for hot reload functionality, initialized when hot reload is enabled
        self._palette_source: str | None = None
        self._watcher: QFileSystemWatcher | None = None
        self._reload_timer: QTimer | None = None

        # Signal for hot reload notifications
        self.signals = ThemeSignals(self)

    def get_theme(self) -> str:
        """Return the current stylesheet content."""
        return self._theme

    def get_palette(self) -> QPalette:
        """Return the configured palette."""
        return self._palette

    def reload(self) -> None:
        """Reload the stylesheet and palette configuration from disk."""
        try:
            updated_theme = (self.theme_dir / _THEME_FILE_NAME).read_text()
        except OSError:
            logger.debug("Theme reload deferred because a theme file is temporarily unavailable")
            if self._watcher is not None and self._reload_timer is not None:
                self._reload_timer.start()
            return

        try:
            updated_palette_source = (self.theme_dir / _PALETTE_FILE_NAME).read_text()
            palette_module = self._load_palette_module(self.theme_dir / _PALETTE_FILE_NAME, updated_palette_source)
            updated_palette = palette_module.build_application_palette()
        except Exception:
            logger.exception("Theme reload deferred because the palette configuration is unavailable or invalid")
            if self._watcher is not None and self._reload_timer is not None:
                self._reload_timer.start()
            return

        self._refresh_watched_paths()

        updated = False

        if updated_palette_source != self._palette_source:
            self._palette_source = updated_palette_source
            self._palette = updated_palette
            updated = True

        if updated_theme != self._theme:
            self._theme = updated_theme
            updated = True

        if updated:
            logger.info("Reloaded GUI theme files from %s", self.theme_dir)
            self.signals.themes_reloaded.emit()

    def enable_hot_reload(self) -> None:
        """Watch the theme directory and reload stylesheets when files change."""
        if self._palette_source is None:
            self._palette_source = (self.theme_dir / _PALETTE_FILE_NAME).read_text()

        if self._watcher is None:
            self._watcher = QFileSystemWatcher(self)
            self._watcher.fileChanged.connect(self._schedule_reload)
            self._watcher.directoryChanged.connect(self._schedule_reload)

        if self._reload_timer is None:
            self._reload_timer = QTimer(self)
            self._reload_timer.setSingleShot(True)
            self._reload_timer.setInterval(100)
            self._reload_timer.timeout.connect(self.reload)

        self._refresh_watched_paths()
        logger.info("Enabled GUI theme hot reload for development")

    def _schedule_reload(self, _path: str) -> None:
        if self._reload_timer is None:
            return

        self._reload_timer.start()

    def _refresh_watched_paths(self) -> None:
        if self._watcher is None:
            return

        file_paths = [str(self.theme_dir / _THEME_FILE_NAME), str(self.theme_dir / _PALETTE_FILE_NAME)]
        watched_files = set(self._watcher.files())
        watched_directories = set(self._watcher.directories())
        expected_directory = str(self.theme_dir)

        missing_files = [path for path in file_paths if path not in watched_files and Path(path).exists()]
        stale_files = [path for path in watched_files if path not in file_paths]
        missing_directory = expected_directory not in watched_directories and self.theme_dir.exists()
        stale_directories = [path for path in watched_directories if path != expected_directory]

        if stale_files:
            self._watcher.removePaths(stale_files)
        if stale_directories:
            self._watcher.removePaths(stale_directories)
        if missing_files:
            self._watcher.addPaths(missing_files)
        if missing_directory:
            self._watcher.addPath(expected_directory)

    @staticmethod
    def _load_palette_module(module_path: Path, source: str) -> ModuleType:
        module_name = "pairinteraction_gui_theme_palette_runtime"
        spec = importlib.util.spec_from_file_location(module_name, module_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Could not load palette module from {module_path}")

        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        try:
            spec.loader.exec_module(module)
            return module
        finally:
            sys.modules.pop(module_name, None)


theme_manager = ThemeManager()
