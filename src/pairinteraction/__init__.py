# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

# ruff: noqa: E402
from __future__ import annotations

import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable


def _setup_dynamic_libraries() -> None:  # noqa: C901, PLR0915
    import platform
    import sys
    from importlib.metadata import PackageNotFoundError, files, version
    from pathlib import Path
    from warnings import warn

    from pairinteraction._info import Info

    # ---------------------------------------------------------------------------------------
    # Helper functions
    # ---------------------------------------------------------------------------------------
    def is_package_installed(package: str) -> bool:
        """Check whether a package is installed."""
        try:
            version(package)
        except PackageNotFoundError:
            return False
        else:
            return True

    def is_running_under_pyinstaller() -> bool:
        return getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS")

    def find_library_file_in_package(substring: str, package: str) -> Path | None:
        """Find the library file if it got installed together with the package."""
        if (package_files := files(package)) is not None:
            for p in package_files:
                if substring in p.stem:
                    return Path(p.locate()).resolve()
        return None

    def find_library_file_in_directory(substring: str, directory: Path = Path(__file__).parent) -> Path | None:
        """Find the library file if it is in the specified directory."""
        for p in directory.glob("*"):
            if substring in p.stem:
                return p.resolve()
        return None

    def load_candidate(candidate: Path, loader: Callable[[str], object]) -> None:
        try:
            loader(str(candidate))
        except Exception as e:
            warn(f"Unable to load {candidate.name}: {e}", RuntimeWarning, stacklevel=2)

    add_dll_directory: Callable[[str], object] | None = getattr(os, "add_dll_directory", None)

    # ---------------------------------------------------------------------------------------
    # Load shared libraries
    # ---------------------------------------------------------------------------------------

    def fix_ssl() -> None:
        """Fix SSL library loading issues under Windows."""
        assert add_dll_directory is not None

        # If PairInteraction was installed, the SSL library might have been installed together with it
        if (
            is_package_installed("pairinteraction")
            and find_library_file_in_package("ssl", "pairinteraction") is not None
        ):
            return

        # If PairInteraction is running under PyInstaller, the SSL library might have been bundled with the executable
        if is_running_under_pyinstaller() and find_library_file_in_directory("ssl", Path(__file__).parent) is not None:
            return

        # Else, PairInteraction is probably running in development mode and we add a bunch of directories to
        # the DLL search path to avoid loading issues
        possible_dirs = [
            Path.cwd(),
            # look in cwd.parents
            Path.cwd().parent,
            Path.cwd().parent.parent,
            Path.cwd().parent.parent.parent,
            # and __file__.parents[2] (for editable installs)
            Path(__file__).resolve().parent.parent.parent,
            Path(__file__).resolve().parent.parent.parent.parent,
        ]
        possible_paths = [d / "vcpkg_installed" / "x64-windows" / "bin" for d in possible_dirs]
        for path in possible_paths:
            if path.is_dir():
                add_dll_directory(str(path))

    def load_mkl(system: str) -> None:
        import ctypes
        from functools import partial

        if not is_package_installed("mkl"):
            raise RuntimeError("The 'mkl' library is not installed.")

        path = find_library_file_in_package("mkl_core", "mkl")
        if path is None:
            raise RuntimeError("The 'mkl_core' library could not be found.")

        mkl_lib_dir = path.parent

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
            tbb_lib_file: Path | None = None
            if is_package_installed("tbb"):
                tbb_lib_file = find_library_file_in_package("tbb", "tbb")
            if tbb_lib_file is None and is_package_installed("pairinteraction"):
                tbb_lib_file = find_library_file_in_package("tbb", "pairinteraction")
            if tbb_lib_file is None:
                raise RuntimeError("The 'tbb' library could not be found.")
            load_candidate(tbb_lib_file, partial(ctypes.CDLL, mode=os.RTLD_LAZY | os.RTLD_GLOBAL))

            for lib in mkl_lib_file_names:
                candidate = mkl_lib_dir / f"lib{lib}.so.2"
                load_candidate(candidate, partial(ctypes.CDLL, mode=os.RTLD_LAZY | os.RTLD_GLOBAL))

        elif system == "Windows":
            assert add_dll_directory is not None

            # Modify the dll search path
            add_dll_directory(str(mkl_lib_dir))

        else:
            warn(f"Cannot load MKL libraries on unsupported system {system}.", RuntimeWarning, stacklevel=2)

    system = platform.system()
    if system == "Windows":
        fix_ssl()
    if Info.with_mkl:
        load_mkl(system)


_setup_dynamic_libraries()


def _setup_ca_bundle() -> None:
    import certifi

    from pairinteraction._backend import set_ca_bundle_path

    set_ca_bundle_path(certifi.where())


_setup_ca_bundle()


# ---------------------------------------------------------------------------------------
# Configure PairInteraction for running tests with a local database if requested
# ---------------------------------------------------------------------------------------
def _setup_test_mode(download_missing: bool = False, database_dir: str | None = None) -> None:
    from pathlib import Path

    from pairinteraction.database import Database

    if database_dir is None:
        possible_dirs = [
            Path.cwd(),
            # look in cwd.parents
            Path.cwd().parent,
            Path.cwd().parent.parent,
            Path.cwd().parent.parent.parent,
            # and __file__.parents[2] (for editable installs)
            Path(__file__).resolve().parent.parent.parent,
        ]
        possible_paths = [d / "data" / "database" for d in possible_dirs]

        for path in possible_paths:
            if any(path.rglob("wigner.parquet")):
                database_dir = str(path)
                break

        if database_dir is None:
            raise FileNotFoundError("Could not find database directory")

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
del _setup_ca_bundle

# ---------------------------------------------------------------------------------------
# Import pairinteraction
# ---------------------------------------------------------------------------------------
from pairinteraction import (
    green_tensor,
    perturbative,
    real,
    visualization,
)
from pairinteraction._backend import (
    VERSION_MAJOR as _VERSION_MAJOR,
    VERSION_MINOR as _VERSION_MINOR,
    VERSION_PATCH as _VERSION_PATCH,
    run_unit_tests,
)
from pairinteraction.basis import BasisAtom, BasisPair
from pairinteraction.custom_logging import configure_logging
from pairinteraction.database import Database, print_database_info
from pairinteraction.diagonalization import diagonalize
from pairinteraction.ket import KetAtom, KetPair
from pairinteraction.perturbative import C3, C6, EffectiveSystemPair
from pairinteraction.state import StateAtom, StatePair
from pairinteraction.system import SystemAtom, SystemPair
from pairinteraction.units import ureg

__all__ = [
    "C3",
    "C6",
    "BasisAtom",
    "BasisPair",
    "Database",
    "EffectiveSystemPair",
    "KetAtom",
    "KetPair",
    "StateAtom",
    "StatePair",
    "SystemAtom",
    "SystemPair",
    "configure_logging",
    "diagonalize",
    "green_tensor",
    "perturbative",
    "print_database_info",
    "real",
    "run_unit_tests",
    "ureg",
    "visualization",
]

__version__ = f"{_VERSION_MAJOR}.{_VERSION_MINOR}.{_VERSION_PATCH}"


# ---------------------------------------------------------------------------------------
# Clean up namespace
# ---------------------------------------------------------------------------------------
del _VERSION_MAJOR, _VERSION_MINOR, _VERSION_PATCH
del _setup_dynamic_libraries  # don't delete _setup_test_mode, since it is used in tests/conftest.py
