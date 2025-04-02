# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import sys

from pairinteraction_gui.app import Application
from pairinteraction_gui.main_window import MainWindow


def main() -> int:
    """Run the PairInteraction GUI application.

    Returns:
        int: Application exit code

    """
    app = Application(sys.argv)
    app.setApplicationName("PairInteraction")

    app.allow_ctrl_c()

    window = MainWindow()
    window.show()

    return sys.exit(app.exec())


if __name__ == "__main__":
    main()
