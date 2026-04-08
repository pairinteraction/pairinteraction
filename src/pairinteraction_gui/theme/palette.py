# SPDX-FileCopyrightText: 2026 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from PySide6.QtGui import QColor, QGuiApplication, QPalette


def build_application_palette() -> QPalette:
    palette = QGuiApplication.palette()

    # Main application surfaces.
    # - Window: main window background, toolbox content areas, toolbox container,
    #   and informational outcome labels
    # - Button: sidebar-adjacent controls including plot toolbar buttons,
    #   toolbox tabs, and regular push buttons
    palette.setColor(QPalette.ColorRole.Window, QColor("#f9fbfe"))
    palette.setColor(QPalette.ColorRole.Button, QColor("#dbe4ef"))

    # Default text on light surfaces.
    # - Text: toolbox page content and informational outcome labels
    # - WindowText: main window text and plot toolbar button text
    # - ButtonText: toolbox tab labels and regular push button text
    # - HighlightedText: text on selected toolbox tabs, pressed plot toolbar
    #   buttons, pressed push buttons, and selected export menu items
    palette.setColor(QPalette.ColorRole.Text, QColor("#111828"))
    palette.setColor(QPalette.ColorRole.WindowText, QColor("#111828"))
    palette.setColor(QPalette.ColorRole.ButtonText, QColor("#111828"))
    palette.setColor(QPalette.ColorRole.HighlightedText, QColor("#111828"))

    # Accent colors for interactive controls.
    # - Highlight: selected toolbox tabs and selected export menu items
    # - Accent: pressed/checked plot toolbar buttons and pressed push buttons
    palette.setColor(QPalette.ColorRole.Highlight, QColor("#cad9f7"))
    palette.setColor(QPalette.ColorRole.Accent, QColor("#cad9f7"))

    # Text intended for dark chrome and emphasized controls.
    # Used for the status bar plus hovered and checked sidebar tool buttons.
    palette.setColor(QPalette.ColorRole.BrightText, QColor("#ffffff"))

    # Neutral support colors for borders, dark panels, and hover states.
    # - Midlight: hover fill for buttons and borders for the toolbox,
    #   info/error labels
    # - Mid: sidebar tool button text and borders in their resting state
    # - Dark: status bar fill, sidebar toolbar background, and menu border
    palette.setColor(QPalette.ColorRole.Midlight, QColor("#cdd5e0"))
    palette.setColor(QPalette.ColorRole.Mid, QColor("#8a91a0"))
    palette.setColor(QPalette.ColorRole.Dark, QColor("#28354e"))

    return palette
