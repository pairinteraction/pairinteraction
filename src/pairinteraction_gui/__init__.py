# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later


import os
import sys

from PySide6.QtCore import QLoggingCategory, QTimer

from pairinteraction_gui.app import Application, SplashScreen

__all__ = ["main"]


def main(*, enable_theme_hot_reload: bool = False) -> int:
    """Run the PairInteraction GUI application.

    Returns:
        int: Application exit code

    """
    rules = os.environ.get("QT_LOGGING_RULES", "")
    if "qt.gui.icc.warning=" not in rules:
        QLoggingCategory.setFilterRules(f"{rules}\nqt.gui.icc.warning=false".strip())

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
        app.processEvents()
        QTimer.singleShot(0, lambda: splash.finish(window))

    # 50 ms lets the compositor flush the first frame, avoiding painting issues with the splash screen on some platforms
    QTimer.singleShot(50, _start)

    return app.exec()
