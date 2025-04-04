# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from functools import wraps
from multiprocessing import Process, SimpleQueue
from typing import Any, Callable

from PySide6.QtCore import QObject, QThread, Signal

from pairinteraction_gui.app import Application

logger = logging.getLogger(__name__)


ALL_PROCESSES: set[Process] = set()
ALL_THREADS: set[QThread] = set()


class WorkerSignals(QObject):
    """Signals to be used by the Worker class."""

    finished = Signal(bool)
    error = Signal(Exception)
    result = Signal(object)


class Worker(QThread):
    """Simple worker class to run a function in a separate thread."""

    def __init__(self, fn: Callable[..., Any], *args: Any, **kwargs: Any) -> None:
        super().__init__(Application.instance())

        ALL_THREADS.add(self)

        self.fn = fn
        self.args = args
        self.kwargs = kwargs

        self.signals = WorkerSignals()
        self.finished.connect(self.finish_up)

    def run(self) -> None:
        """Initialise the runner function with passed args, kwargs."""
        logger.debug("Starting thread %s", self)
        success = False
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
        ALL_THREADS.remove(self)


def run_in_other_process(func: Callable[..., Any]) -> Callable[..., Any]:
    @wraps(func)
    def wrapper_func(*args: Any, **kwargs: Any) -> Any:
        def mp_func(queue: SimpleQueue, *args: Any, **kwargs: Any) -> None:
            result = func(*args, **kwargs)
            queue.put(result)

        queue = SimpleQueue()
        process = Process(target=mp_func, args=(queue, *args), kwargs=kwargs, daemon=True)
        ALL_PROCESSES.add(process)
        process.start()
        logger.debug("Starting process %s", process.pid)
        result = queue.get()
        process.join()
        logger.debug("Closing process %s", process.pid)
        process.close()
        ALL_PROCESSES.remove(process)
        return result

    return wrapper_func
