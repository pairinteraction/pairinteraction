# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later


import multiprocessing
import sys

from pairinteraction_gui.app import Application
from pairinteraction_gui.main_window import MainWindow

__all__ = ["main"]


def main() -> int:
    """Run the PairInteraction GUI application.

    Returns:
        int: Application exit code

    """
    # Multithreading together with "fork" is not supported.
    # Furthermore, "spawn" will become the default in Python 3.14 on all platforms,
    # see also https://docs.python.org/3/library/multiprocessing.html#contexts-and-start-methods
    # Thus, we already now set the start method to "spawn" for all platforms.
    multiprocessing.set_start_method("spawn")

    app = Application(sys.argv)
    app.setApplicationName("PairInteraction")

    app.allow_ctrl_c()

    window = MainWindow()
    window.show()

    return app.exec()
