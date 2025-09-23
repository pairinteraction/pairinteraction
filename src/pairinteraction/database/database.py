# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, ClassVar

from pairinteraction import _backend
from pairinteraction.custom_logging import _flush_pending_logs

if TYPE_CHECKING:
    import os
    from pathlib import Path


logger = logging.getLogger(__name__)


class Database:
    """Class for handling the databases for the PairInteraction package.

    The Databases are used to store the atomic states, their energies, and their matrix elements to other states.
    The databases are stored in the user's cache directory by default, but can be stored in a different directory.
    When running PairInteraction for the first time, the databases have to be downloaded from the internet
    (e.g. by explicitly passing `download_missing=True` to the constructor).
    Once the databases are downloaded, the user usually does not have to interact with the Database class directly.

    """

    _global_database: ClassVar[Database | None] = None

    def __init__(
        self,
        download_missing: bool = False,
        use_cache: bool = True,
        database_dir: str | os.PathLike[str] = "",
    ) -> None:
        """Create a new database instance with the given parameters.

        Args:
            download_missing: Whether to download missing databases if needed. Default False.
            use_cache: Whether to load the Wigner 3j symbols table into memory. Default True.
            database_dir: The directory where the databases are stored.
                Default "", i.e. use the default directory (the user's cache directory).

        """
        self._cpp = _backend.Database(download_missing, use_cache, database_dir)
        _flush_pending_logs()  # call it manually because constructors of classes from nanobind cannot be decorated

    @classmethod
    def _from_cpp_object(cls, cpp_obj: _backend.Database) -> Database:
        """Create a Database instance from a C++ Database object.

        This is used internally to convert C++ objects returned by the C++ API to Python objects.
        """
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    @classmethod
    def get_global_database(cls) -> Database:
        """Return the global database instance if it was initialized, otherwise None."""
        return cls._global_database  # type: ignore [return-value]

    @classmethod
    def initialize_global_database(
        cls,
        download_missing: bool = False,
        use_cache: bool = True,
        database_dir: str | os.PathLike[str] = "",
    ) -> None:
        """Initialize the global database with the given parameters.

        The arguments are the same as for the constructor of this class.
        """
        db = cls(download_missing, use_cache, database_dir)
        if cls._global_database is None:
            cls._global_database = db
        elif (
            cls._global_database.download_missing == db.download_missing
            and cls._global_database.use_cache == db.use_cache
            and cls._global_database.database_dir == db.database_dir
        ):
            pass  # already initialized with the same parameters
        else:
            raise ValueError(
                "Global database was already initialized with different parameters. "
                "The global database is automatically initialized when needed. "
                "If you explicitly want to initialize the global database, do this at the beginning of your script."
            )

    @property
    def download_missing(self) -> bool:
        """Whether to download missing databases if needed."""
        return self._cpp.get_download_missing()

    @property
    def use_cache(self) -> bool:
        """Whether to load the Wigner 3j symbols table into memory."""
        return self._cpp.get_use_cache()

    @property
    def database_dir(self) -> Path:
        """The directory where the databases are stored."""
        return self._cpp.get_database_dir()
