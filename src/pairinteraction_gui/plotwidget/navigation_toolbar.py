# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import Optional

from matplotlib.backends.backend_qt import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from PySide6.QtWidgets import QWidget


class CustomNavigationToolbar(NavigationToolbar):
    """Custom navigation toolbar for matplotlib figures.

    See Also:
    https://stackoverflow.com/questions/12695678/how-to-modify-the-navigation-toolbar-easily-in-a-matplotlib-figure-window/15549675#15549675

    """

    toolitems = (
        ("Zoom", "Zoom to rectangle\nx/y fixes axis", "zoom_to_rect", "zoom"),
        ("Pan", "Left button pans, Right button zooms\nx/y fixes axis, CTRL fixes aspect", "move", "pan"),
        ("Home", "Reset original view", "home", "home"),
    )

    def __init__(self, canvas: FigureCanvasQTAgg, parent: Optional[QWidget] = None) -> None:
        """Initialize the custom navigation toolbar."""
        super().__init__(canvas, parent, coordinates=False)
