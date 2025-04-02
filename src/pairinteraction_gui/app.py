# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import os
import signal
from types import FrameType
from typing import Optional

from PySide6.QtCore import QCoreApplication, QObject, QSocketNotifier, QTimer, Signal
from PySide6.QtWidgets import QApplication


class MainSignals(QObject):
    """Signals for the application.

    We store an instance of this signal class in the Application instance, see app.py.
    So to access these signals (from anywhere in the application), you can use
    `Application.instance().signals`.
    """

    ask_download_database = Signal(str)


class Application(QApplication):
    """Add some global signals to the QApplication."""

    signals = MainSignals()

    @classmethod
    def instance(cls) -> "Application":  # type: ignore  # overwrite type hints
        """Return the current instance of the application."""
        return super().instance()  # type: ignore [return-value]

    def allow_ctrl_c(self) -> None:
        # Create a pipe to communicate between the signal handler and the Qt event loop
        pipe_r, pipe_w = os.pipe()

        def signal_handler(signal: int, frame: Optional[FrameType]) -> None:
            os.write(pipe_w, b"x")  # Write a single byte to the pipe

        signal.signal(signal.SIGINT, signal_handler)

        def handle_signal() -> None:
            os.read(pipe_r, 1)  # Read the byte from the pipe to clear it
            print("\nCtrl+C detected. Shutting down gracefully...")
            QCoreApplication.quit()

        sn = QSocketNotifier(pipe_r, QSocketNotifier.Type.Read, parent=self)
        sn.activated.connect(handle_signal)

        # Create a timer to ensure the event loop processes events regularly
        # This makes Ctrl+C work even when the application is idle
        timer = QTimer(self)
        timer.timeout.connect(lambda: None)  # Do nothing, just wake up the event loop
        timer.start(200)
