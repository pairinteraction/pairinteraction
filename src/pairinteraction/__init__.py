# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

# ruff: noqa: E402

import os
from typing import Optional


def _setup_dynamic_libaries() -> None:  # noqa: C901, PLR0915
    import os
    import platform
    from importlib.metadata import PackageNotFoundError, files
    from pathlib import Path
    from typing import Callable
    from warnings import warn

    from pairinteraction._info import Info

    # ---------------------------------------------------------------------------------------
    # Helper functions
    # ---------------------------------------------------------------------------------------
    def get_library_file(package: str, substring: str) -> Path:
        try:
            package_files = files(package)
            if package_files is None:
                raise RuntimeError(f"Installation database records of '{package}' are missing.")
            file_path = Path(next(p for p in package_files if substring in p.stem).locate())
            return file_path.resolve()
        except PackageNotFoundError as err:
            # Also look in the current directory for the library (used by pyinstaller)
            for p in Path(__file__).parent.glob("*"):
                if substring in p.stem:
                    return p.resolve()
            raise RuntimeError(f"The '{package}' package could not be found.") from err
        except StopIteration as err:
            raise RuntimeError(f"The '{substring}' library could not be found.") from err

    def load_candidate(candidate: Path, loader: Callable[..., None]) -> None:
        try:
            loader(str(candidate))
        except Exception as e:
            warn(f"Unable to load {candidate.name}: {e}", RuntimeWarning, stacklevel=2)

    # ---------------------------------------------------------------------------------------
    # Load shared libraries
    # ---------------------------------------------------------------------------------------
    def fix_ssl() -> None:
        try:
            get_library_file("pairinteraction", "ssl")
        except RuntimeError as err:
            if "library could not be found" not in str(err):
                raise
        else:
            return

        vcpkg_dirs = [
            Path.cwd() / "vcpkg_installed/x64-windows/bin",
            Path.cwd().parent / "vcpkg_installed/x64-windows/bin",
            Path.cwd().parent.parent / "vcpkg_installed/x64-windows/bin",
        ]
        if directory := next((d for d in vcpkg_dirs if d.is_dir()), None):
            os.add_dll_directory(str(directory))  # type: ignore [attr-defined]

    def load_mkl(system: str) -> None:
        import ctypes
        import platform
        import sys
        from functools import partial

        mkl_lib_dir = get_library_file("mkl", "mkl_core").parent

        mkl_lib_file_names = [
            "mkl_core",  # must be loaded first
            "mkl_tbb_thread",  # must be loaded second
            "mkl_avx2",
            "mkl_avx512",
            "mkl_def",
            "mkl_intel_lp64",
            "mkl_mc3",
            "mkl_rt",
            "mkl_vml_avx2",
            "mkl_vml_avx512",
            "mkl_vml_cmpt",
            "mkl_vml_def",
            "mkl_vml_mc3",
            "mkl_sequential",  # needed for pytest with cpp coverage
        ]

        if system == "Linux":
            # Under linux, the libraries must always be loaded manually in the address space
            try:
                tbb_lib_file = get_library_file("tbb", "tbb")
            except RuntimeError:
                tbb_lib_file = get_library_file("pairinteraction", "tbb")
            load_candidate(tbb_lib_file, partial(ctypes.CDLL, mode=os.RTLD_LAZY | os.RTLD_GLOBAL))  # type: ignore [arg-type]

            for lib in mkl_lib_file_names:
                candidate = mkl_lib_dir / f"lib{lib}.so.2"
                load_candidate(candidate, partial(ctypes.CDLL, mode=os.RTLD_LAZY | os.RTLD_GLOBAL))  # type: ignore [arg-type]

        elif system == "Windows":
            # Modify the dll search path
            os.add_dll_directory(str(mkl_lib_dir))  # type: ignore [attr-defined]

            is_conda_cpython = platform.python_implementation() == "CPython" and (
                hasattr(ctypes.pythonapi, "Anaconda_GetVersion") or "packaged by conda-forge" in sys.version
            )

            if sys.version_info[:2] >= (3, 10) or not is_conda_cpython:
                return

            # If conda python <= 3.9 is used, the libraries must be loaded manually in the address space
            # https://github.com/adang1345/delvewheel/blob/c37a82f0f66dd73e0169ff637f7c0ba5b33032c6/delvewheel/_wheel_repair.py#L56-L77
            from ctypes import WinDLL, wintypes  # type: ignore [attr-defined]

            kernel32 = WinDLL("kernel32", use_last_error=True)
            kernel32.LoadLibraryExW.restype = wintypes.HMODULE
            kernel32.LoadLibraryExW.argtypes = wintypes.LPCWSTR, wintypes.HANDLE, wintypes.DWORD
            for lib in mkl_lib_file_names:
                candidate = mkl_lib_dir / f"{lib}.2.dll"
                load_candidate(candidate, partial(kernel32.LoadLibraryExW, None, 8))

        else:
            warn(f"Cannot load MKL libraries on unsupported system {system}.", RuntimeWarning, stacklevel=2)

    system = platform.system()
    if system == "Windows":
        fix_ssl()
    if Info.with_mkl:
        load_mkl(system)


_setup_dynamic_libaries()


# ---------------------------------------------------------------------------------------
# Configure pairinteraction for running tests with a local database if requested
# ---------------------------------------------------------------------------------------
def _setup_test_mode(download_missing: bool = False, database_dir: Optional[str] = None) -> None:
    from pathlib import Path

    from pairinteraction._wrapped import Database

    pairinteraction_dir = Path(__file__).parent
    cwd = Path.cwd()

    if database_dir is None:
        possible_dirs = [
            pairinteraction_dir.parent.parent / "data" / "database",  # if local editable install
            cwd / "data" / "database",  # normal pytest mode
            cwd.parent / "data" / "database",  # for pytest jupyter notebooks
            cwd.parent.parent / "data" / "database",  # for pytest jupyter notebooks
            cwd.parent.parent.parent / "data" / "database",  # for pytest jupyter notebooks
        ]
        try:
            database_dir = next(str(d) for d in possible_dirs if any(d.rglob("wigner.parquet")))
        except StopIteration:
            raise FileNotFoundError("Could not find database directory") from None

    Database.initialize_global_database(download_missing, True, database_dir)


if os.getenv("PAIRINTERACTION_TEST_MODE", "0") == "1":
    download_missing = bool(int(os.getenv("PAIRINTERACTION_TEST_DOWNLOAD_MISSING", "0")))
    database_dir = os.getenv("PAIRINTERACTION_TEST_DATABASE_DIR", None)
    _setup_test_mode(download_missing, database_dir)


# ---------------------------------------------------------------------------------------
# Decorate all functions in _backend with a decorator that flushes pending logs
# ---------------------------------------------------------------------------------------
def _setup_logging() -> None:
    from pairinteraction import _backend
    from pairinteraction.custom_logging import decorate_module_with_flush_logs

    decorate_module_with_flush_logs(_backend)


_setup_logging()
del _setup_logging

# ---------------------------------------------------------------------------------------
# Import pairinteraction
# ---------------------------------------------------------------------------------------
from pairinteraction import (
    complex,  # noqa: A004
    perturbative,
    real,
    visualization,
)
from pairinteraction._backend import (
    VERSION_MAJOR as _VERSION_MAJOR,
    VERSION_MINOR as _VERSION_MINOR,
    VERSION_PATCH as _VERSION_PATCH,
)
from pairinteraction.custom_logging import configure_logging
from pairinteraction.module_tests import run_module_tests
from pairinteraction.units import ureg

__all__ = [
    "complex",
    "configure_logging",
    "perturbative",
    "real",
    "run_module_tests",
    "ureg",
    "visualization",
]

__version__ = f"{_VERSION_MAJOR}.{_VERSION_MINOR}.{_VERSION_PATCH}"


# ---------------------------------------------------------------------------------------
# Clean up namespace
# ---------------------------------------------------------------------------------------
del _VERSION_MAJOR, _VERSION_MINOR, _VERSION_PATCH
del _setup_dynamic_libaries  # don't delete _setup_test_mode, since it is used in tests/conftest.py
