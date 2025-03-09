import datetime
import inspect
import logging
import re
from typing import Any

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
    created = datetime.datetime.strptime(fields["timestamp"], "%Y-%m-%d %H:%M:%S.%f").timestamp()
    record = logging.LogRecord(
        name=logger.name,
        level=level,
        pathname=fields["filename"],
        lineno=int(fields["lineno"]),
        msg=fields["message"],
        args=(),
        exc_info=None,
        created=created,
        thread=int(fields["thread"]),
    )
    if level >= logger.getEffectiveLevel():
        logger.handle(record)


def _flush_pending_logs() -> None:
    for entry in get_pending_logs():
        _log_cpp_backend_record(entry.level, entry.message)


def _flush_logs_after(func: callable) -> callable:
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
