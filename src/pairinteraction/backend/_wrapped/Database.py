from typing import TYPE_CHECKING, ClassVar, Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.OperatorType import OperatorType, get_cpp_operator_type

if TYPE_CHECKING:
    import os

    from scipy.sparse import csr_matrix

    from pairinteraction.backend._wrapped.basis.BasisAtom import BasisAtom

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
    def get_global_instance(cls) -> "Database":
        if not hasattr(cls, "GlobalDatabase"):
            cls.create_global_instance()
        return cls.GlobalDatabase

    @classmethod
    def create_global_instance(
        cls,
        download_missing: bool = False,
        wigner_in_memory: bool = True,
        database_dir: Union[str, "os.PathLike[str]"] = "",
    ) -> None:
        if hasattr(cls, "GlobalDatabase"):
            raise ValueError("Global instance already exists")
        # FIXME since CPPDatabase.get_global_instance is currently broken, we use a normal Database as global instance
        cls.GlobalDatabase = cls(download_missing, wigner_in_memory, database_dir)

    def get_matrix_elements(
        self, basis_ket: "BasisAtom", basis_bra: "BasisAtom", operator: OperatorType, q: int
    ) -> "csr_matrix":
        cpp_operator_type = get_cpp_operator_type(operator)
        return self._cpp.get_matrix_elements(basis_ket._cpp, basis_bra._cpp, cpp_operator_type, q)  # type: ignore
