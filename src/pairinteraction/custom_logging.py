# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import datetime
import inspect
import logging
import re
from typing import Any, Callable, ClassVar

from colorama import Fore, Style, just_fix_windows_console

from pairinteraction._backend import get_pending_logs


def _extract_cpp_backend_log_fields(message: str) -> dict[str, str]:
    pattern = (
        r"^\[(?P<timestamp>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d+)\s+(?P<thread>\d+)\]\s*"
        r".*"
        r"\[(?P<filename>[^:\]]+):(?P<lineno>\d+)\]\s*"
        r"(?P<message>.*)$"
    )
    match = re.match(pattern, message.strip(), re.DOTALL)
    if not match:
        raise RuntimeError(f"Could not parse log message: {message}")
    return match.groupdict()


def _log_cpp_backend_record(level: int, message: str) -> None:
    logger = logging.getLogger("cpp")
    fields = _extract_cpp_backend_log_fields(message)
    record = logging.LogRecord(
        name=logger.name,
        level=level,
        pathname=fields["filename"],
        lineno=int(fields["lineno"]),
        msg=fields["message"],
        args=(),
        exc_info=None,
    )
    record.created = datetime.datetime.strptime(fields["timestamp"], "%Y-%m-%d %H:%M:%S.%f").timestamp()
    record.thread = int(fields["thread"])
    if level >= logger.getEffectiveLevel():
        logger.handle(record)


def _flush_pending_logs() -> None:
    for entry in get_pending_logs():
        _log_cpp_backend_record(entry.level, entry.message)


def _flush_logs_after(func: Callable[..., Any]) -> Callable[..., Any]:
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        result = func(*args, **kwargs)
        _flush_pending_logs()
        return result

    return wrapper


def decorate_module_with_flush_logs(module: object) -> None:
    for name, obj in vars(module).items():
        if inspect.isclass(obj):
            for attr_name, attr in vars(obj).items():
                if callable(attr) and not attr_name.startswith("__"):
                    setattr(obj, attr_name, _flush_logs_after(attr))
        elif callable(obj) and not name.startswith("__"):
            setattr(module, name, _flush_logs_after(obj))


def configure_logging(
    level_str: str = "WARNING",
    fmt: str = ("[%(asctime)s.%(msecs)03d] [%(levelname)s] [%(filename)s:%(lineno)d] %(message)s"),
) -> None:
    """Configure colorfully formatted logging."""

    class ColoredFormatter(logging.Formatter):
        COLORS: ClassVar = {
            "DEBUG": Fore.BLUE,
            "INFO": Fore.GREEN,
            "WARNING": Fore.YELLOW,
            "ERROR": Fore.RED,
            "CRITICAL": Fore.RED + Style.BRIGHT,
        }

        def format(self, record: logging.LogRecord) -> str:
            original_levelname = record.levelname
            record.levelname = f"{self.COLORS[record.levelname]}{record.levelname}{Style.RESET_ALL}"
            formatted = super().format(record)
            record.levelname = original_levelname
            return formatted

    just_fix_windows_console()

    handler = logging.StreamHandler()
    handler.setFormatter(ColoredFormatter(fmt, datefmt="%H:%M:%S"))

    root_logger = logging.getLogger()
    if root_logger.hasHandlers():
        root_logger.handlers.clear()
    root_logger.setLevel(getattr(logging, level_str.upper()))
    root_logger.addHandler(handler)
