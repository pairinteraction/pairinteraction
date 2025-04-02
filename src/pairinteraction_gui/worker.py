# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from functools import wraps
from typing import Any, Callable

from PySide6.QtCore import QObject, QRunnable, QThreadPool, Signal, Slot

THREADPOOL = QThreadPool()
THREADPOOL.setMaxThreadCount(2)


class Worker(QRunnable):
    """Simple worker class to run a function in a separate thread."""

    def __init__(self, fn: Callable[..., Any], *args: Any, **kwargs: Any) -> None:
        super().__init__()

        self.fn = fn
        self.args = args
        self.kwargs = kwargs

        self.signals = WorkerSignals()

    @Slot()
    def run(self) -> None:
        """Initialise the runner function with passed args, kwargs."""
        try:
            result = self.fn(*self.args, **self.kwargs)
        except Exception:
            import sys
            import traceback

            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit()


class WorkerSignals(QObject):
    """Signals to be used by the Worker class."""

    finished = Signal()
    error = Signal(tuple)
    result = Signal(object)
    progress = Signal(int)


def run_worker(fn: Callable[..., Any], *args: Any, **kwargs: Any) -> None:
    """Run a function in a separate thread."""
    worker = Worker(fn, *args, **kwargs)
    THREADPOOL.start(worker)


def threaded(func: Callable[..., Any]) -> Callable[..., Any]:
    @wraps(func)
    def wrapper_func(*args: Any, **kwargs: Any) -> Any:
        worker = Worker(func, *args, **kwargs)
        THREADPOOL.start(worker)
        return worker.signals.result

    return wrapper_func


# def threaded_factory(func_finished: Optional[Callable[..., Any]] = None) -> Callable[..., Any]:
#     """Decorator to run the function in a separate thread."""

#     def threaded(func: Callable[..., Any]) -> Callable[..., Any]:
#         @wraps(func)
#         def wrapper_func(*args: Any, **kwargs: Any) -> Any:
#             worker = Worker(func, *args, **kwargs)
#             if func_finished:
#                 worker.signals.finished.connect(func_finished)
#             THREADPOOL.start(worker)
#             return worker.signals.result

#         return wrapper_func

#     return threaded


"""
def example_usage():
    worker = Worker(execute_this_fn, "foo", "bar", baz=42)
    worker.signals.result.connect(print_output)
    worker.signals.finished.connect(thread_complete)
    worker.signals.progress.connect(progress_fn)
    THREADPOOL.start(worker)
"""
