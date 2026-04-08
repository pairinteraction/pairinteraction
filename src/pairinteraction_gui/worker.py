# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, ClassVar

from PySide6.QtCore import QObject, QSize, Qt, QThread, Signal
from PySide6.QtGui import QMovie
from PySide6.QtWidgets import QApplication, QLabel

if TYPE_CHECKING:
    from PySide6.QtWidgets import QWidget


logger = logging.getLogger(__name__)


class WorkerSignals(QObject):
    """Signals to be used by the Worker class."""

    started = Signal()
    finished = Signal(bool)
    error = Signal(Exception)
    result = Signal(object)


class MultiThreadWorker(QThread):
    """Simple worker class to run a function in a separate thread.

    Example:
    worker = Worker(my_function, arg1, arg2, kwarg1=value1)
    worker.signals.result.connect(process_result)
    worker.start()

    """

    all_threads: ClassVar[set[MultiThreadWorker]] = set()

    def __init__(self, fn: Callable[..., Any], *args: Any, **kwargs: Any) -> None:
        super().__init__(QApplication.instance())

        self.all_threads.add(self)

        self.fn = fn
        self.args = args
        self.kwargs = kwargs

        self.signals = WorkerSignals()
        self.finished.connect(self.finish_up)

    def enable_busy_indicator(self, widget: QWidget) -> None:
        """Run a loading gif while the worker is running."""
        self.busy_label = QLabel(widget)
        gif_path = Path(__file__).parent / "images" / "loading.gif"
        self.busy_movie = QMovie(str(gif_path))
        self.busy_movie.setScaledSize(QSize(100, 100))  # Make the gif larger
        self.busy_label.setMovie(self.busy_movie)
        self.busy_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.busy_label.setGeometry((widget.width() - 100) // 2, (widget.height() - 100) // 2, 100, 100)

        self.signals.started.connect(self._start_gif)
        self.signals.finished.connect(self._stop_gif)

    def _start_gif(self) -> None:
        self.busy_label.show()
        self.busy_movie.start()

    def _stop_gif(self, _success: bool = True) -> None:
        if hasattr(self, "busy_movie"):
            self.busy_movie.stop()
        if hasattr(self, "busy_label"):
            self.busy_label.hide()

    def run(self) -> None:
        """Initialise the runner function with passed args, kwargs."""
        logger.debug("Run on thread %s", self)
        success = False
        self.signals.started.emit()
        try:
            result = self.fn(*self.args, **self.kwargs)
            success = True
        except Exception as err:
            self.signals.error.emit(err)
        else:
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit(success)

    def finish_up(self) -> None:
        """Perform any final cleanup or actions before the thread exits."""
        logger.debug("Finishing up thread %s", self)
        self.all_threads.discard(self)

    @classmethod
    def terminate_all(cls) -> None:
        """Terminate all threads started by the application."""
        # Shallow copy to avoid error if the set is modified during the loop,
        # e.g. if the thread is finished and removes itself from the list
        all_threads = list(cls.all_threads)
        for thread in all_threads:
            if thread.isRunning():
                logger.debug("Terminating thread %s.", thread)
                thread.terminate()
                thread.signals.finished.emit(False)
                thread.wait()

        cls.all_threads.clear()
        logger.debug("All threads terminated.")
