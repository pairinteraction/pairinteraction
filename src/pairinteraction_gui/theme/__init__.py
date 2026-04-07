# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from pathlib import Path
from types import ModuleType
from typing import Protocol, cast

from PySide6.QtCore import QFileSystemWatcher, QObject, QTimer, Signal
from PySide6.QtGui import QPalette

logger = logging.getLogger(__name__)

_THEME_FILE_NAME = "theme.qss"
_PALETTE_FILE_NAME = "palette.py"


class ThemeSignals(QObject):
    """Signals emitted when theme files change."""

    themes_reloaded = Signal()


class PaletteModule(Protocol):
    def build_application_palette(self) -> QPalette: ...


class ThemeManager(QObject):
    """Load theme files and optionally watch them for development reloads."""

    def __init__(self, theme_dir: Path | None = None) -> None:
        super().__init__()
        self.signals = ThemeSignals(self)
        self.theme_dir = theme_dir or Path(__file__).parent
        self._theme = ""
        self._palette_source = ""
        self._palette = QPalette()
        self._watcher: QFileSystemWatcher | None = None
        self._reload_timer = QTimer(self)
        self._reload_timer.setSingleShot(True)
        self._reload_timer.setInterval(75)
        self._reload_timer.timeout.connect(self.reload)
        self.reload(emit_signal=False)

    def get_theme(self) -> str:
        """Return the current stylesheet content."""
        return self._theme

    def get_palette(self) -> QPalette:
        """Return a copy of the configured palette."""
        return QPalette(self._palette)

    def reload(self, *, emit_signal: bool = True) -> None:
        """Reload the stylesheet and palette configuration from disk."""
        try:
            updated_theme = (self.theme_dir / _THEME_FILE_NAME).read_text()
            updated_palette_source = (self.theme_dir / _PALETTE_FILE_NAME).read_text()
        except OSError:
            logger.debug("Theme reload deferred because a theme file is temporarily unavailable")
            if self._watcher is not None:
                self._reload_timer.start()
            return

        try:
            palette_module = self._load_palette_module(self.theme_dir / _PALETTE_FILE_NAME, updated_palette_source)
            updated_palette = palette_module.build_application_palette()
        except Exception:
            logger.exception("Theme reload deferred because the palette configuration is invalid")
            if self._watcher is not None:
                self._reload_timer.start()
            return

        if updated_theme == self._theme and updated_palette_source == self._palette_source:
            self._refresh_watched_paths()
            return

        self._theme = updated_theme
        self._palette_source = updated_palette_source
        self._palette = updated_palette
        self._refresh_watched_paths()
        logger.info("Reloaded GUI theme files from %s", self.theme_dir)
        if emit_signal:
            self.signals.themes_reloaded.emit()

    def enable_hot_reload(self) -> None:
        """Watch the theme directory and reload stylesheets when files change."""
        if self._watcher is None:
            self._watcher = QFileSystemWatcher(self)
            self._watcher.fileChanged.connect(self._schedule_reload)
            self._watcher.directoryChanged.connect(self._schedule_reload)

        self._refresh_watched_paths()
        logger.info("Enabled GUI theme hot reload for development")

    def set_theme_dir(self, theme_dir: Path) -> None:
        """Update the theme directory. Intended for tests and development tooling."""
        self.theme_dir = theme_dir
        self.reload()

    def _schedule_reload(self, _path: str) -> None:
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
    def _load_palette_module(module_path: Path, source: str) -> PaletteModule:
        module = ModuleType("pairinteraction_gui_theme_palette_runtime")
        module.__file__ = str(module_path)
        exec(compile(source, str(module_path), "exec"), module.__dict__)  # noqa: S102
        return cast("PaletteModule", module)


theme_manager = ThemeManager()
