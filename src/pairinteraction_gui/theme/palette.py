# SPDX-FileCopyrightText: 2026 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from PySide6.QtGui import QColor, QPalette


def build_application_palette() -> QPalette:
    palette = QPalette()

    # Main application surfaces.
    # - Window: main window background
    # - Base: content areas such as toolbox pages and info labels
    # - Button: push buttons
    palette.setColor(QPalette.ColorRole.Window, QColor("#ffffff"))
    palette.setColor(QPalette.ColorRole.Base, QColor("#f8fafc"))
    palette.setColor(QPalette.ColorRole.Button, QColor("#e4e8ee"))

    # Default text on light surfaces.
    # Used for regular content text, main window text, and button labels.
    palette.setColor(QPalette.ColorRole.Text, QColor("#111828"))
    palette.setColor(QPalette.ColorRole.WindowText, QColor("#111828"))
    palette.setColor(QPalette.ColorRole.ButtonText, QColor("#111828"))

    # Accent color for selected and active elements.
    # Used for checked sidebar buttons, selected toolbox tabs,
    # and checked buttons in the plot navigation toolbar.
    palette.setColor(QPalette.ColorRole.Highlight, QColor("#cad9f7"))
    palette.setColor(QPalette.ColorRole.HighlightedText, QColor("#111828"))

    # Text intended for dark or accent-colored backgrounds.
    # Used for the status bar, sidebar toolbar buttons, and toolbox tabs.
    palette.setColor(QPalette.ColorRole.BrightText, QColor("#ffffff"))

    # Neutral support colors for borders, dark panels, and hover states.
    # - Light: bright border/accent for toolbox tabs
    # - Midlight: soft borders, including info label outlines
    # - Mid: hover background for interactive controls
    # - Dark: dark chrome such as the status bar, sidebar toolbar, and toolbox tabs
    palette.setColor(QPalette.ColorRole.Light, QColor("#ffffff"))
    palette.setColor(QPalette.ColorRole.Midlight, QColor("#cdd5e0"))
    palette.setColor(QPalette.ColorRole.Mid, QColor("#8a91a0"))
    palette.setColor(QPalette.ColorRole.Dark, QColor("#28354e"))

    return palette
