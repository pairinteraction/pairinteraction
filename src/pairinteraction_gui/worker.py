# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
import os
from functools import wraps
from multiprocessing.pool import Pool
from pathlib import Path
from threading import Thread
from typing import TYPE_CHECKING, Any, Callable, ClassVar, Optional, TypeVar

from PySide6.QtCore import QObject, QSize, Qt, QThread, Signal
from PySide6.QtGui import QMovie
from PySide6.QtWidgets import QApplication, QLabel, QWidget

if TYPE_CHECKING:
    from typing_extensions import ParamSpec

    P = ParamSpec("P")
    R = TypeVar("R")

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

    all_threads: ClassVar[set["MultiThreadWorker"]] = set()

    def __init__(self, fn: Callable[..., Any], *args: Any, **kwargs: Any) -> None:
        super().__init__(QApplication.instance())

        self.all_threads.add(self)

        self.fn = fn
        self.args = args
        self.kwargs = kwargs

        self.signals = WorkerSignals()
        self.finished.connect(self.finish_up)

    def enable_busy_indicator(self, widget: "QWidget") -> None:
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


class MultiProcessWorker:
    _mp_functions_dict: ClassVar[dict[str, Callable[..., Any]]] = {}
    _pool: ClassVar[Optional[Pool]] = None
    _async_worker: ClassVar[Optional[Thread]] = None

    def __init__(self, fn_name: str, *args: Any, **kwargs: Any) -> None:
        if fn_name not in self._mp_functions_dict:
            raise ValueError(f"Function {fn_name} is not registered.")

        self.fn_name = fn_name
        self.args = args
        self.kwargs = kwargs

    @classmethod
    def create_pool(cls, n_processes: int = 1) -> None:
        """Create a pool of processes."""
        if cls._pool is not None or cls._async_worker is not None:
            raise RuntimeError(
                "create_pool already called. Use terminate_all(create_new_pool=True) to restart the pool."
            )

        cls._async_worker = Thread(target=cls._create_pool, args=(n_processes,))
        cls._async_worker.start()

    @classmethod
    def _create_pool(cls, n_processes: int) -> None:
        """Create a pool of processes."""
        cls._pool = Pool(n_processes)
        cls._pool.apply(cls._dummy_function)  # Call the pool once, to make the next call faster
        cls._async_worker = None
        logger.debug("Pool created successfully.")

    @staticmethod
    def _dummy_function() -> None:
        """Do nothing.

        Dummy function to run after creating the pool asynchronously.
        """
        return

    @classmethod
    def register(cls, func: Callable[..., Any], name: Optional[str] = None) -> None:
        name = name if name is not None else func.__name__
        if name in cls._mp_functions_dict:
            raise ValueError(f"Function {name} is already registered.")
        cls._mp_functions_dict[name] = func

    def start(self) -> Any:
        async_worker = self._async_worker
        if async_worker is not None:
            logger.debug("Waiting for creating_pool to finish.")
            async_worker.join()
            logger.debug("creating_pool finished.")

        if self._pool is None:
            raise RuntimeError("Pool is not created. Call create_pool() first.")

        logger.debug("Starting pool.apply")
        result = self._pool.apply(self.run)
        logger.debug("Finished pool.apply")
        return result

    def run(self) -> Any:
        logger.debug("Run on process %s", os.getpid())
        func = self._mp_functions_dict[self.fn_name]
        return func(*self.args, **self.kwargs)

    @classmethod
    def terminate_all(cls, create_new_pool: bool) -> None:
        """Terminate all processes."""
        if cls._async_worker is not None:
            cls._async_worker.join()

        if cls._pool is None:
            return

        cls._pool.terminate()
        cls._pool = None
        logger.debug("Process pool terminated.")

        if create_new_pool:
            cls.create_pool()


def run_in_other_process(func: Callable["P", "R"]) -> Callable["P", "R"]:
    MultiProcessWorker.register(func)

    @wraps(func)
    def wrapper_func(*args: "P.args", **kwargs: "P.kwargs") -> "R":
        return MultiProcessWorker(func.__name__, *args, **kwargs).start()  # type: ignore [no-any-return]

    return wrapper_func
