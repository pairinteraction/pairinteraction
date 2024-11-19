from typing import TYPE_CHECKING, ClassVar, Union

from pairinteraction.backend import _backend

if TYPE_CHECKING:
    import os

CPPDatabase = _backend.Database


class Database:
    GlobalDatabase: ClassVar["Database"]

    def __init__(
        self,
        download_missing: bool = False,
        wigner_in_memory: bool = True,
        database_dir: Union[str, "os.PathLike[str]"] = "",
    ) -> None:
        self._cpp = CPPDatabase(download_missing, wigner_in_memory, database_dir)

    @classmethod
    def _from_cpp_object(cls, cpp_obj: CPPDatabase) -> "Database":
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    @classmethod
    def get_global_instance(
        cls,
        download_missing: bool = False,
        wigner_in_memory: bool = True,
        database_dir: Union[str, "os.PathLike[str]"] = "",
    ) -> "Database":
        if not hasattr(cls, "GlobalDatabase"):
            cls.GlobalDatabase = cls._from_cpp_object(
                CPPDatabase.get_global_instance(download_missing, wigner_in_memory, database_dir)  # type: ignore [reportUnkownMemberType]
            )
        return cls.GlobalDatabase
