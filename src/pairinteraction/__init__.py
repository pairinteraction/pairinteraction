import os
import platform
from typing import Union


# ---------------------------------------------------------------------------------------
# Load shared libraries
# ---------------------------------------------------------------------------------------
def _fix_ssl() -> None:
    from importlib.metadata import files
    from pathlib import Path

    if next((p for p in files("pairinteraction") if "ssl" in p.stem), None):
        return

    vcpkg_dirs = [
        Path.cwd() / "vcpkg_installed/x64-windows/bin",
        Path.cwd().parent / "vcpkg_installed/x64-windows/bin",
        Path.cwd().parent.parent / "vcpkg_installed/x64-windows/bin",
    ]
    if directory := next((d for d in vcpkg_dirs if d.is_dir()), None):
        os.add_dll_directory(str(directory))


def _load_mkl(system: str) -> None:
    import ctypes
    import platform
    import sys
    from functools import partial
    from importlib.metadata import PackageNotFoundError, files
    from pathlib import Path

    def get_library_file(package: str, substring: str) -> Path:
        try:
            return next(p for p in files(package) if substring in p.stem).locate().resolve()
        except (PackageNotFoundError, StopIteration):
            raise RuntimeError(f"The '{substring}' library could not be found.") from None

    def load_candidate(candidate: os.PathLike, loader: callable) -> None:
        if candidate.exists():
            loader(str(candidate))

    tbb_lib_file = get_library_file("pairinteraction", "tbb")
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
    ]

    if system == "Linux":
        # Under linux, the libraries must always be loaded manually in the address space
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
        os.add_dll_directory(str(mkl_lib_dir))
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.LoadLibraryExW.restype = ctypes.wintypes.HMODULE
        kernel32.LoadLibraryExW.argtypes = ctypes.wintypes.LPCWSTR, ctypes.wintypes.HANDLE, ctypes.wintypes.DWORD
        for lib in mkl_lib_file_names:
            candidate = mkl_lib_dir / f"{lib}.2.dll"
            load_candidate(candidate, partial(kernel32.LoadLibraryExW, None, 8))


system = platform.system()
if system == "Windows":
    _fix_ssl()
    _load_mkl(system)
elif system == "Linux":
    _load_mkl(system)


# ---------------------------------------------------------------------------------------
# Import pairinteraction
# ---------------------------------------------------------------------------------------
from pairinteraction import backend  # noqa: E402
from pairinteraction.backend._backend import VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH  # noqa: E402
from pairinteraction.module_tests import run_module_tests  # noqa: E402
from pairinteraction.units import ureg  # noqa: E402

__all__ = [
    "backend",
    "run_module_tests",
    "ureg",
]

__version__ = f"{VERSION_MAJOR}.{VERSION_MINOR}.{VERSION_PATCH}"


# ---------------------------------------------------------------------------------------
# Configure pairinteraction for running tests with a local database if requested
# ---------------------------------------------------------------------------------------
def setup_test_mode(download_missing: bool = False, database_dir: Union[str, "os.PathLike[str]", None] = None) -> None:
    from pathlib import Path

    from pairinteraction.backend._wrapped import Database

    PAIRINTERACTION_DIR = Path(__file__).parent
    CWD = Path.cwd()

    if database_dir is None:
        possible_dirs = [
            PAIRINTERACTION_DIR.parent.parent / "data" / "database",  # if local editable install
            CWD / "data" / "database",  # normal pytest mode
            CWD.parent / "data" / "database",  # for pytest jupyter notebooks
            CWD.parent.parent / "data" / "database",  # for pytest jupyter notebooks
            CWD.parent.parent.parent / "data" / "database",  # for pytest jupyter notebooks
        ]
        if not (database_dir := next((d for d in possible_dirs if any(d.glob("wigner*.parquet"))), None)):
            raise FileNotFoundError("Could not find database directory")
    else:
        database_dir = Path(database_dir)

    Database.initialize_global_database(download_missing, True, database_dir)


if os.getenv("PAIRINTERACTION_TEST_MODE", "0") == "1":
    setup_test_mode()

# Remove symbols imported for internal use.
# Note that we do not delete setup_test_mode, because it is also used in tests/conftest.py.
del platform, _fix_ssl, _load_mkl
