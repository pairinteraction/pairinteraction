# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from pathlib import Path

__all__ = ["label_error_theme", "label_theme", "main_theme"]


cwd = Path(__file__).parent
main_theme = (cwd / "main.qss").read_text()
label_theme = (cwd / "label.qss").read_text()
label_error_theme = (cwd / "label_error.qss").read_text()
