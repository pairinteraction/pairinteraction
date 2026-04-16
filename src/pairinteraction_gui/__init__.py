# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later


import sys

from PySide6.QtCore import QTimer

from pairinteraction_gui.app import Application, SplashScreen

__all__ = ["main"]


def main(*, enable_theme_hot_reload: bool = False) -> int:
    """Run the PairInteraction GUI application.

    Returns:
        int: Application exit code

    """
    app = Application(sys.argv)
    app.setApplicationName("PairInteraction")
    app.allow_ctrl_c()

    splash = SplashScreen()
    splash.show()
    app.processEvents()

    def _start() -> None:
        from pairinteraction_gui.main_window import MainWindow

        window = MainWindow(enable_theme_hot_reload=enable_theme_hot_reload)
        window.show()
        splash.finish(window)

    # 50 ms lets the compositor flush the first frame, avoiding painting issues with the splash screen on some platforms
    QTimer.singleShot(50, _start)

    return app.exec()
