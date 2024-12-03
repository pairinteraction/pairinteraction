from typing import TYPE_CHECKING, ClassVar, Union, overload

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.OperatorType import OperatorType, get_cpp_operator_type
from pairinteraction.units import QuantitySparse

if TYPE_CHECKING:
    import os

    from pint.facets.plain import PlainQuantity
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

    @overload
    def get_matrix_elements(
        self, basis_ket: "BasisAtom", basis_bra: "BasisAtom", operator: OperatorType, q: int
    ) -> "PlainQuantity[csr_matrix]": ...  # type: ignore [reportInvalidTypeArguments]

    @overload
    def get_matrix_elements(
        self, basis_ket: "BasisAtom", basis_bra: "BasisAtom", operator: OperatorType, q: int, unit: str
    ) -> "csr_matrix": ...

    def get_matrix_elements(
        self, basis_ket: "BasisAtom", basis_bra: "BasisAtom", operator: OperatorType, q: int, unit: str = "pint"
    ):
        cpp_operator_type = get_cpp_operator_type(operator)
        matrix_elements_au = self._cpp.get_matrix_elements(basis_ket._cpp, basis_bra._cpp, cpp_operator_type, q)  # type: ignore
        matrix_elements = QuantitySparse.from_base(matrix_elements_au, operator)
        return matrix_elements.to_unit(unit)
