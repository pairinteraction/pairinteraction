import os
import platform
from pathlib import Path


# ---------------------------------------------------------------------------------------
# Load the MKL library from the mkl pypi package
# ---------------------------------------------------------------------------------------
def _load_mkl() -> None:
    import ctypes
    import sys
    from functools import partial
    from importlib.metadata import PackageNotFoundError, files

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

    lib_file_names = [
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

    if platform.system() == "Linux":
        # Under linux, the libraries must always be loaded manually in the address space
        load_candidate(tbb_lib_file, partial(ctypes.CDLL, mode=os.RTLD_LAZY | os.RTLD_GLOBAL))
        for lib in lib_file_names:
            candidate = mkl_lib_dir / f"lib{lib}.so.2"
            load_candidate(candidate, partial(ctypes.CDLL, mode=os.RTLD_LAZY | os.RTLD_GLOBAL))

    else:
        # Modify the dll search path
        os.add_dll_directory(str(mkl_lib_dir))

        vcpkg_dirs = [
            Path.cwd() / "vcpkg_installed/x64-windows/bin",
            Path.cwd().parent / "vcpkg_installed/x64-windows/bin",
            Path.cwd().parent.parent / "vcpkg_installed/x64-windows/bin",
        ]
        if directory := next((d for d in vcpkg_dirs if d.is_dir()), None):
            os.add_dll_directory(str(directory))

        is_conda_cpython = platform.python_implementation() == "CPython" and (
            hasattr(ctypes.pythonapi, "Anaconda_GetVersion") or "packaged by conda-forge" in sys.version
        )

        if sys.version_info[:2] >= (3, 10) or not is_conda_cpython:
            return

        # If conda python <= 3.9 is used, the libraries must be loaded manually in the address space
        os.add_dll_directory(str(mkl_lib_dir))
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.LoadLibraryExW.restype = ctypes.wintypes.HMODULE
        kernel32.LoadLibraryExW.argtypes = ctypes.wintypes.LPCWSTR, ctypes.wintypes.HANDLE, ctypes.wintypes.DWORD
        for lib in lib_file_names:
            candidate = mkl_lib_dir / f"{lib}.2.dll"
            load_candidate(candidate, partial(kernel32.LoadLibraryExW, None, 8))


if platform.system() in ["Linux", "Windows"]:
    _load_mkl()
del _load_mkl

# ---------------------------------------------------------------------------------------
# Import pairinteraction
# ---------------------------------------------------------------------------------------
from pairinteraction import backend  # noqa: E402
from pairinteraction.backend._backend import VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH  # noqa: E402
from pairinteraction.module_tests import run_module_tests  # noqa: E402
from pairinteraction.setup_test_mode import setup_test_mode  # noqa: E402
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
if os.getenv("PAIRINTERACTION_TEST_MODE", "0") == "1":
    setup_test_mode()
