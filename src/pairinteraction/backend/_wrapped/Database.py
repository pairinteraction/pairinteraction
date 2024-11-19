from typing import TYPE_CHECKING, Union

import pairinteraction.backend._backend as _backend

if TYPE_CHECKING:
    import os

CPPDatabase = _backend.Database

GlobalDatabase = None


class Database:
    def __init__(
        self,
        download_missing: bool = False,
        wigner_in_memory: bool = True,
        database_dir: Union[str, "os.PathLike"] = "",
    ) -> None:
        self._cpp = CPPDatabase(download_missing, wigner_in_memory, database_dir)

    @classmethod
    def _from_cpp_object(cls, cpp_obj: CPPDatabase):
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    @staticmethod
    def get_global_instance(
        download_missing: bool = False, wigner_in_memory: bool = True, database_dir: Union[str, "os.PathLike"] = ""
    ) -> "Database":
        global GlobalDatabase
        if GlobalDatabase is None:
            GlobalDatabase = Database._from_cpp_object(
                CPPDatabase.get_global_instance(download_missing, wigner_in_memory, database_dir)
            )
        return GlobalDatabase
