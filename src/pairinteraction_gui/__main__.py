import logging
import sys

from pairinteraction_gui.app import Application
from pairinteraction_gui.main_window import MainWindow


def main() -> int:
    """Run the PairInteraction GUI application.

    Returns:
        int: Application exit code

    """
    # in case we run the app from the terminal, add some logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    app = Application(sys.argv)
    app.setApplicationName("PairInteraction")

    app.allow_ctrl_c()

    window = MainWindow()
    window.show()

    return sys.exit(app.exec())


if __name__ == "__main__":
    main()
