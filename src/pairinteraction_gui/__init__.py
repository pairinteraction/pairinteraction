# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later


import sys

from pairinteraction_gui.app import Application
from pairinteraction_gui.main_window import MainWindow

__all__ = ["main"]


def main(*, enable_theme_hot_reload: bool = False) -> int:
    """Run the PairInteraction GUI application.

    Returns:
        int: Application exit code

    """
    app = Application(sys.argv)
    app.setApplicationName("PairInteraction")

    app.allow_ctrl_c()

    window = MainWindow(enable_theme_hot_reload=enable_theme_hot_reload)
    window.show()

    return app.exec()
