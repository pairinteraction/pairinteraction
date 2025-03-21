# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from functools import wraps
from multiprocessing import Process, SimpleQueue
from typing import TYPE_CHECKING, Any, Callable, TypeVar

from PySide6.QtCore import QObject, QThread, Signal

from pairinteraction_gui.app import Application

if TYPE_CHECKING:
    from typing_extensions import ParamSpec

    P = ParamSpec("P")
    R = TypeVar("R")

logger = logging.getLogger(__name__)


class WorkerSignals(QObject):
    """Signals to be used by the Worker class."""

    finished = Signal(bool)
    error = Signal(Exception)
    result = Signal(object)


class Worker(QThread):
    """Simple worker class to run a function in a separate thread."""

    def __init__(self, fn: Callable[..., Any], *args: Any, **kwargs: Any) -> None:
        super().__init__(Application.instance())

        Application.all_threads.add(self)

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
        Application.all_threads.remove(self)


def run_in_other_process(func: Callable["P", "R"]) -> Callable["P", "R"]:
    @wraps(func)
    def wrapper_func(*args: "P.args", **kwargs: "P.kwargs") -> "R":
        def mp_func(queue: "SimpleQueue[R]", *args: "P.args", **kwargs: "P.kwargs") -> None:
            result = func(*args, **kwargs)
            queue.put(result)

        queue: SimpleQueue[R] = SimpleQueue()
        process = Process(target=mp_func, args=(queue, *args), kwargs=kwargs, daemon=True)
        Application.all_processes.add(process)
        process.start()
        logger.debug("Starting process %s", process.pid)
        result = queue.get()
        process.join()
        logger.debug("Closing process %s", process.pid)
        process.close()
        Application.all_processes.remove(process)
        return result

    return wrapper_func
