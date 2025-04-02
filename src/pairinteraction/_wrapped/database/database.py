# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from typing import TYPE_CHECKING, ClassVar, Optional, Union

from pairinteraction import _backend

if TYPE_CHECKING:
    import os


logger = logging.getLogger(__name__)

CPPDatabase = _backend.Database


class Database:
    """Class for a handling the databases for the pairinteraction package.

    The Databases are used to store the atomic states, their energies, and their matrix elements to other states.
    The databases are stored in the user's cache directory by default, but can be stored in a different directory.
    When running pairinteraction for the first time, the databases have to be downloaded from the internet
    (e.g. by explicitly passing `download_missing=True` to the constructor).
    Once the databases are downloaded, the user usually does not have to interact with the Database class directly.

    """

    _global_database: ClassVar[Optional["Database"]] = None

    def __init__(
        self,
        download_missing: bool = False,
        use_cache: bool = True,
        database_dir: Union[str, "os.PathLike[str]"] = "",
    ) -> None:
        """Create a new database instance with the given parameters.

        Args:
            download_missing: Whether to download missing databases if needed. Default False.
            use_cache: Whether to load the Wigner 3j symbols table into memory. Default True.
            database_dir: The directory where the databases are stored.
                Default "", i.e. use the default directory (the user's cache directory).

        """
        self._cpp = CPPDatabase(download_missing, use_cache, database_dir)
        self.download_missing = download_missing
        self.use_cache = use_cache
        self.database_dir = database_dir

    @classmethod
    def get_global_database(cls) -> "Database":
        """Return the global database instance if it was initialized, otherwise None."""
        return cls._global_database  # type: ignore [return-value]

    @classmethod
    def initialize_global_database(
        cls,
        download_missing: bool = False,
        use_cache: bool = True,
        database_dir: Union[str, "os.PathLike[str]"] = "",
    ) -> None:
        """Initialize the global database with the given parameters.

        The arguments are the same as for the constructor of this class.
        """
        if cls._global_database is None:
            cls._global_database = cls(download_missing, use_cache, database_dir)
        elif (
            cls._global_database.download_missing == download_missing
            and cls._global_database.use_cache == use_cache
            and cls._global_database.database_dir == database_dir
        ):
            pass  # already initialized with the same parameters
        else:
            raise ValueError(
                "Global database was already initialized with different parameters. "
                "The global database is automatically initialized when needed. "
                "If you explicitly want to initialize the global database, do this at the beginning of your script."
            )
