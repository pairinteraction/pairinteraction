def _setup_dynamic_libaries() -> None:  # noqa: C901
    import os
    import platform
    from importlib.metadata import PackageNotFoundError, files
    from pathlib import Path
    from typing import Callable
    from warnings import warn

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
        except (PackageNotFoundError, StopIteration):
            raise RuntimeError(f"The '{substring}' library could not be found.") from None

    def load_candidate(candidate: os.PathLike[str], loader: Callable) -> None:
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
            return
        except RuntimeError as err:
            if "library could not be found" not in str(err):
                raise err

        vcpkg_dirs = [
            Path.cwd() / "vcpkg_installed/x64-windows/bin",
            Path.cwd().parent / "vcpkg_installed/x64-windows/bin",
            Path.cwd().parent.parent / "vcpkg_installed/x64-windows/bin",
        ]
        if directory := next((d for d in vcpkg_dirs if d.is_dir()), None):
            os.add_dll_directory(str(directory))

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
            load_candidate(tbb_lib_file, partial(ctypes.CDLL, mode=os.RTLD_LAZY | os.RTLD_GLOBAL))

            for lib in mkl_lib_file_names:
                candidate = mkl_lib_dir / f"lib{lib}.so.2"
                load_candidate(candidate, partial(ctypes.CDLL, mode=os.RTLD_LAZY | os.RTLD_GLOBAL))

        elif system == "Windows":
            # Modify the dll search path
            os.add_dll_directory(str(mkl_lib_dir))

            is_conda_cpython = platform.python_implementation() == "CPython" and (
                hasattr(ctypes.pythonapi, "Anaconda_GetVersion") or "packaged by conda-forge" in sys.version
            )

            if sys.version_info[:2] >= (3, 10) or not is_conda_cpython:
                return

            # If conda python <= 3.9 is used, the libraries must be loaded manually in the address space
            # https://github.com/adang1345/delvewheel/blob/c37a82f0f66dd73e0169ff637f7c0ba5b33032c6/delvewheel/_wheel_repair.py#L56-L77
            import ctypes.wintypes as wintypes

            kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
            kernel32.LoadLibraryExW.restype = wintypes.HMODULE
            kernel32.LoadLibraryExW.argtypes = wintypes.LPCWSTR, wintypes.HANDLE, wintypes.DWORD
            for lib in mkl_lib_file_names:
                candidate = mkl_lib_dir / f"{lib}.2.dll"
                load_candidate(candidate, partial(kernel32.LoadLibraryExW, None, 8))

    system = platform.system()
    if system == "Windows":
        fix_ssl()
        load_mkl(system)
    elif system == "Linux":
        load_mkl(system)
    elif system == "Darwin":
        pass
    else:
        warn(f"Cannot load MKL libraries on unknown system {system}.", RuntimeWarning, stacklevel=2)


_setup_dynamic_libaries()


# ---------------------------------------------------------------------------------------
# Configure pairinteraction for running tests with a local database if requested
# ---------------------------------------------------------------------------------------
def _setup_test_mode(download_missing: bool = False, database_dir: str = "") -> None:
    from pathlib import Path

    from pairinteraction._wrapped import Database

    PAIRINTERACTION_DIR = Path(__file__).parent
    CWD = Path.cwd()

    if database_dir == "":
        possible_dirs = [
            PAIRINTERACTION_DIR.parent.parent / "data" / "database",  # if local editable install
            CWD / "data" / "database",  # normal pytest mode
            CWD.parent / "data" / "database",  # for pytest jupyter notebooks
            CWD.parent.parent / "data" / "database",  # for pytest jupyter notebooks
            CWD.parent.parent.parent / "data" / "database",  # for pytest jupyter notebooks
        ]
        try:
            database_dir = next(str(d) for d in possible_dirs if any(d.glob("wigner*.parquet")))
        except StopIteration:
            raise FileNotFoundError("Could not find database directory") from None

    Database.initialize_global_database(download_missing, True, Path(database_dir))


import os as _os  # noqa: E402

if _os.getenv("PAIRINTERACTION_TEST_MODE", "0") == "1":
    _setup_test_mode()


# ---------------------------------------------------------------------------------------
# Import pairinteraction
# ---------------------------------------------------------------------------------------
from pairinteraction import complex, real  # noqa: E402
from pairinteraction._backend import (  # noqa: E402
    VERSION_MAJOR as _VERSION_MAJOR,
    VERSION_MINOR as _VERSION_MINOR,
    VERSION_PATCH as _VERSION_PATCH,
)
from pairinteraction.module_tests import run_module_tests  # noqa: E402
from pairinteraction.units import ureg  # noqa: E402

__all__ = [
    "complex",
    "real",
    "run_module_tests",
    "ureg",
]

__version__ = f"{_VERSION_MAJOR}.{_VERSION_MINOR}.{_VERSION_PATCH}"


# ---------------------------------------------------------------------------------------
# Clean up namespace
# ---------------------------------------------------------------------------------------
del _VERSION_MAJOR, _VERSION_MINOR, _VERSION_PATCH
del _os
del _setup_dynamic_libaries  #  don't delete _setup_test_mode, since it is used in tests/conftest.py
