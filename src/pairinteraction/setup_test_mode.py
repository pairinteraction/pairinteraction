import os
from pathlib import Path
from typing import Union

from pairinteraction.backend._wrapped import Database

PAIRINTERACTION_DIR = Path(__file__).parent
CWD = Path.cwd()


def setup_test_mode(download_missing: bool = False, database_dir: Union[str, "os.PathLike[str]", None] = None) -> None:
    if database_dir is None:
        possible_dirs = [
            PAIRINTERACTION_DIR.parent.parent / "data" / "database",  # if local editable install
            CWD / "data" / "database",  # normal pytest mode
            CWD.parent / "data" / "database",  # for pytest jupyter notebooks
            CWD.parent.parent / "data" / "database",  # for pytest jupyter notebooks
        ]

        for database_dir in possible_dirs:
            if any(database_dir.glob("wigner*.parquet")):
                break
        else:
            raise FileNotFoundError("Could not find database directory")
    else:
        database_dir = Path(database_dir)

    Database.initialize_global_database(download_missing, True, database_dir)
