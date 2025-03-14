import logging
import os
import signal
import sys
from types import FrameType
from typing import Optional

from PySide6.QtCore import QCoreApplication, QSocketNotifier, QTimer
from PySide6.QtWidgets import QApplication

from pairinteraction_gui.main_window import MainWindow


def main() -> int:
    """Run the PairInteraction GUI application.

    Returns:
        int: Application exit code

    """
    # in case we run the app from the terminal, add some logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    app = QApplication(sys.argv)
    app.setApplicationName("PairInteraction")

    allow_ctrl_c(app)

    window = MainWindow()
    window.show()

    return sys.exit(app.exec())


def allow_ctrl_c(app: QApplication) -> None:
    # # Make program killable via Ctrl+C if it is started in a terminal
    # signal.signal(signal.SIGINT, signal.SIG_DFL)  # this kills the program

    # This is a more reliable way to handle Ctrl+C in Qt applications
    # Create a pipe to communicate between the signal handler and the Qt event loop
    pipe_r, pipe_w = os.pipe()

    def signal_handler(signal: int, frame: Optional[FrameType]) -> None:
        os.write(pipe_w, b"x")  # Write a single byte to the pipe

    signal.signal(signal.SIGINT, signal_handler)

    def handle_signal() -> None:
        os.read(pipe_r, 1)  # Read the byte from the pipe to clear it
        print("\nCtrl+C detected. Shutting down gracefully...")
        QCoreApplication.quit()

    sn = QSocketNotifier(pipe_r, QSocketNotifier.Type.Read, parent=app)
    sn.activated.connect(handle_signal)

    # Create a timer to ensure the event loop processes events regularly
    # This makes Ctrl+C work even when the application is idle
    timer = QTimer(app)
    timer.timeout.connect(lambda: None)  # Do nothing, just wake up the event loop
    timer.start(200)


if __name__ == "__main__":
    main()
