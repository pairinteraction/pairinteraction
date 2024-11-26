from typing import TYPE_CHECKING, Any, ClassVar, Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.OperatorType import OperatorType, get_cpp_operator_type

if TYPE_CHECKING:
    import os

    import scipy.sparse

    from pairinteraction.backend._wrapped.basis.BasisAtom import BasisAtomBase

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
        self.download_missing = download_missing
        self.wigner_in_memory = wigner_in_memory
        self.database_dir = database_dir

    @classmethod
    def get_global_instance(
        cls,
        download_missing: bool = False,
        wigner_in_memory: bool = True,
        database_dir: Union[str, "os.PathLike[str]"] = "",
    ) -> "Database":
        if not hasattr(cls, "GlobalDatabase"):
            obj = cls.__new__(cls)
            obj._cpp = CPPDatabase.get_global_instance(download_missing, wigner_in_memory, database_dir)  # type: ignore [reportUnkownMemberType]
            obj.download_missing = download_missing
            obj.wigner_in_memory = wigner_in_memory
            obj.database_dir = database_dir
            cls.GlobalDatabase = obj
        if (
            cls.GlobalDatabase.download_missing != download_missing
            or cls.GlobalDatabase.wigner_in_memory != wigner_in_memory
            or cls.GlobalDatabase.database_dir != database_dir
        ):
            raise ValueError("Global instance already exists with different parameters")
        return cls.GlobalDatabase

    def get_matrix_elements(
        self, basis_ket: "BasisAtomBase[Any]", basis_bra: "BasisAtomBase[Any]", operator: OperatorType, q: int
    ) -> "scipy.sparse.csr_matrix":
        cpp_operator_type = get_cpp_operator_type(operator)
        return self._cpp.get_matrix_elements(basis_ket._cpp, basis_bra._cpp, cpp_operator_type, q)  # type: ignore
