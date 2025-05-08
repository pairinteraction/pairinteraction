# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later


def _load_theme_file(file_name: str) -> str:
    """Load a theme file from the current directory."""
    from pathlib import Path

    return (Path(__file__).parent / file_name).read_text()


main_theme = _load_theme_file("main.qss")
label_theme = _load_theme_file("label.qss")
label_error_theme = _load_theme_file("label_error.qss")
plot_widget_theme = _load_theme_file("plot_widget.qss")
