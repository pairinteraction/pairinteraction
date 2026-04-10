# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
import math
from collections.abc import Callable
from typing import TYPE_CHECKING, Any, ClassVar

from PySide6.QtCore import QObject, Qt, QThread, QTimer, Signal
from PySide6.QtGui import QColor, QPainter
from PySide6.QtWidgets import QApplication, QWidget

from pairinteraction import _backend

if TYPE_CHECKING:
    from collections.abc import Callable

logger = logging.getLogger(__name__)


class SpinnerWidget(QWidget):
    """Spinning wheel indicator similar to the OS busy cursor, but larger."""

    def __init__(self, parent: QWidget, size: int = 80, n_dots: int = 12, dot_radius: int = 5) -> None:
        super().__init__(parent)
        self._step = 0
        self._n_dots = n_dots
        self._dot_radius = dot_radius
        self._timer = QTimer(self)
        self._timer.setInterval(80)
        self._timer.timeout.connect(self._advance)
        self.setFixedSize(size, size)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.hide()

    def _advance(self) -> None:
        self._step = (self._step + 1) % self._n_dots
        self.update()

    def start(self) -> None:
        self._timer.start()
        self.show()

    def stop(self) -> None:
        self._timer.stop()
        self.hide()

    def paintEvent(self, event: Any) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        center = self.width() / 2
        orbit = center - self._dot_radius - 4
        for i in range(self._n_dots):
            opacity = ((i - self._step) % self._n_dots) / (self._n_dots - 1)
            painter.setBrush(QColor(128, 128, 128, int(opacity * 220)))
            painter.setPen(Qt.PenStyle.NoPen)
            angle = 2 * math.pi * i / self._n_dots
            x = center + orbit * math.sin(angle)
            y = center - orbit * math.cos(angle)
            painter.drawEllipse(
                int(x - self._dot_radius),
                int(y - self._dot_radius),
                self._dot_radius * 2,
                self._dot_radius * 2,
            )


class WorkerSignals(QObject):
    """Signals to be used by the Worker class."""

    started = Signal()
    finished = Signal(str)
    error = Signal(Exception)
    progress = Signal(str)
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

        self._last_progress_message = ""
        self._progress_timer = QTimer(self)
        self._progress_timer.setInterval(75)
        self._progress_timer.timeout.connect(lambda: self.report_progress(_backend.get_task_status()))
        self.started.connect(self._progress_timer.start)
        self.finished.connect(self._progress_timer.stop)
        self.finished.connect(self.finish_up)

    @classmethod
    def current_worker(cls) -> MultiThreadWorker | None:
        current_thread = QThread.currentThread()
        if isinstance(current_thread, cls):
            return current_thread
        return None

    @classmethod
    def task_checkpoint(cls, progress_message: str | None = None) -> None:
        worker = cls.current_worker()

        if progress_message is not None and worker is not None:
            worker.report_progress(progress_message)

        if worker is not None and worker.isInterruptionRequested():
            raise _backend.TaskAbortedError

        _backend.task_checkpoint(progress_message or "")

    def report_progress(self, message: str) -> None:
        if not message or message == self._last_progress_message:
            return

        self._last_progress_message = message
        self.signals.progress.emit(message)

    def enable_busy_indicator(self, widget: QWidget) -> None:
        """Show a spinning wheel overlay while the worker is running."""
        size = 80
        self.busy_spinner = SpinnerWidget(widget, size=size)
        self.busy_spinner.setGeometry((widget.width() - size) // 2, (widget.height() - size) // 2, size, size)

        self.signals.started.connect(self.busy_spinner.start)
        self.signals.finished.connect(self.busy_spinner.stop)

    def run(self) -> None:
        """Initialise the runner function with passed args, kwargs."""
        logger.debug("Run on thread %s", self)
        status = None
        self.signals.started.emit()
        try:
            result = self.fn(*self.args, **self.kwargs)
            status = "Calculation succeeded"
        except _backend.TaskAbortedError:
            logger.debug("Calculation thread %s aborted.", self)
            status = "Calculation aborted"
        except Exception as err:
            self.signals.error.emit(err)
        else:
            self.signals.result.emit(result)
        finally:
            if status is None:
                status = "Calculation failed"
            _backend.clear_task_abort()
            self.signals.finished.emit(status)

    def request_abort(self) -> None:
        """Request a cooperative abort for the running task."""
        logger.debug("Requesting abort for thread %s.", self)
        self.requestInterruption()
        _backend.request_task_abort()

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
                thread.request_abort()

        for thread in all_threads:
            if thread.isRunning():
                logger.debug("Waiting for thread %s to be aborted.", thread)
                thread.wait()

        cls.all_threads.clear()
        logger.debug("All threads terminated.")
