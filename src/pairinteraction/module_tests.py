import os
from typing import Union

from pairinteraction._backend import run_unit_tests


def run_module_tests(download_missing: bool = False, database_dir: Union[str, "os.PathLike[str]"] = "") -> int:
    wigner_in_memory = True
    return run_unit_tests(download_missing, wigner_in_memory, database_dir)
