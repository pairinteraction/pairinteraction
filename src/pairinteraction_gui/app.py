# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
import os
import signal
from multiprocessing import Process
from types import FrameType
from typing import ClassVar, Optional

from PySide6.QtCore import QObject, QSocketNotifier, QThread, QTimer, Signal
from PySide6.QtWidgets import QApplication

logger = logging.getLogger(__name__)


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
    all_processes: ClassVar[set[Process]] = set()
    all_threads: ClassVar[set[QThread]] = set()

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
            logger.info("Ctrl+C detected in terminal. Shutting down gracefully...")
            self.quit()

        sn = QSocketNotifier(pipe_r, QSocketNotifier.Type.Read, parent=self)
        sn.activated.connect(handle_signal)

        # Create a timer to ensure the event loop processes events regularly
        # This makes Ctrl+C work even when the application is idle
        timer = QTimer(self)
        timer.timeout.connect(lambda: None)  # Do nothing, just wake up the event loop
        timer.start(200)

    @staticmethod
    def quit() -> None:
        """Quit the application."""
        logger.debug("Calling Application.quit().")
        Application.terminate_all_processes()
        Application.terminate_all_threads()
        QApplication.quit()
        logger.debug("Application.quit() done.")

    @staticmethod
    def terminate_all_processes() -> None:
        """Terminate all processes started by the application."""
        for process in Application.all_processes:
            if process.is_alive():
                logger.debug("Terminating process %s.", process.pid)
                process.terminate()
                process.join(timeout=1)

        Application.all_processes.clear()
        logger.debug("All processes terminated.")

    @staticmethod
    def terminate_all_threads() -> None:
        """Terminate all threads started by the application."""
        for thread in Application.all_threads:
            if thread.isRunning():
                logger.debug("Terminating thread %s.", thread)
                thread.terminate()
                thread.wait()

        Application.all_threads.clear()
        logger.debug("All threads terminated.")
