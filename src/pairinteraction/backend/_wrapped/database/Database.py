import logging
from typing import TYPE_CHECKING, ClassVar, Optional, Union, overload

import numpy as np

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.OperatorType import OperatorType, get_cpp_operator_type
from pairinteraction.units import QuantityScalar, QuantitySparse

if TYPE_CHECKING:
    import os

    from pint.facets.plain import PlainQuantity
    from scipy.sparse import csr_matrix

    from pairinteraction.backend._wrapped.basis.BasisAtom import BasisAtom
    from pairinteraction.backend._wrapped.ket.KetAtom import KetAtom
    from pairinteraction.backend._wrapped.system.SystemAtom import SystemAtom

logger = logging.getLogger(__name__)

CPPDatabase = _backend.Database


class Database:
    _global_database: ClassVar[Optional["Database"]] = None

    def __init__(
        self,
        download_missing: bool = False,
        wigner_in_memory: bool = True,
        database_dir: Union[str, "os.PathLike[str]"] = "",
    ) -> None:
        """Create a new database instance with the given parameters.

        Args:
            download_missing: Whether to download missing databases if needed. Default False.
            wigner_in_memory: Whether to load the Wigner 3j symbols table into memory. Default True.
            database_dir: The directory where the databases are stored.
                Default "", i.e. use the default directory (the user's cache directory).

        """
        self._cpp = CPPDatabase(download_missing, wigner_in_memory, database_dir)
        self.download_missing = download_missing
        self.wigner_in_memory = wigner_in_memory
        self.database_dir = database_dir

    @classmethod
    def get_global_database(cls) -> "Database":
        """Return the global database instance if it was initialized, otherwise None."""
        return cls._global_database

    @classmethod
    def initialize_global_database(
        cls,
        download_missing: bool = False,
        wigner_in_memory: bool = True,
        database_dir: Union[str, "os.PathLike[str]"] = "",
    ) -> None:
        """Initialize the global database with the given parameters.

        The arguments are the same as for the constructor of this class.
        """
        if cls._global_database is None:
            cls._global_database = cls(download_missing, wigner_in_memory, database_dir)
        elif (
            cls._global_database.download_missing == download_missing
            and cls._global_database.wigner_in_memory == wigner_in_memory
            and cls._global_database.database_dir == database_dir
        ):
            pass  # already initialized with the same parameters
        else:
            raise ValueError(
                "Global database was already initialized with different parameters. "
                "The global database is automatically initialized when needed. "
                "If you explicitly want to initialize the global database, do this at the beginning of your script."
            )

    @overload
    def get_matrix_elements(
        self, basis_1: "BasisAtom", basis_2: "BasisAtom", operator: OperatorType, q: int
    ) -> "PlainQuantity[csr_matrix]": ...  # type: ignore [reportInvalidTypeArguments]

    @overload
    def get_matrix_elements(
        self, basis_1: "BasisAtom", basis_2: "BasisAtom", operator: OperatorType, q: int, unit: str
    ) -> "csr_matrix": ...

    def get_matrix_elements(
        self, basis_1: "BasisAtom", basis_2: "BasisAtom", operator: OperatorType, q: int, unit: str = "pint"
    ):
        cpp_operator_type = get_cpp_operator_type(operator)
        matrix_elements_au = self._cpp.get_matrix_elements(basis_1._cpp, basis_2._cpp, cpp_operator_type, q)  # type: ignore
        matrix_elements = QuantitySparse.from_base(matrix_elements_au, operator)
        return matrix_elements.to_unit(unit)

    @overload
    def get_matrix_element(
        self, ket_1: "KetAtom", ket_2: "KetAtom", operator: OperatorType, q: int
    ) -> "PlainQuantity[float]": ...  # type: ignore [reportInvalidTypeArguments]

    @overload
    def get_matrix_element(
        self, ket_1: "KetAtom", ket_2: "KetAtom", operator: OperatorType, q: int, *, unit: str
    ) -> "float": ...

    @overload
    def get_matrix_element(
        self, ket_1: "KetAtom", ket_2: "KetAtom", operator: OperatorType, q: int, *, system: "SystemAtom"
    ) -> "PlainQuantity[float]": ...  # type: ignore [reportInvalidTypeArguments]

    @overload
    def get_matrix_element(
        self, ket_1: "KetAtom", ket_2: "KetAtom", operator: OperatorType, q: int, *, system: "SystemAtom", unit: str
    ) -> "float": ...

    def get_matrix_element(
        self,
        ket_1: "KetAtom",
        ket_2: "KetAtom",
        operator: OperatorType,
        q: int,
        *,
        system: Optional["SystemAtom"] = None,
        unit: str = "pint",
    ):
        if system is None:
            BasisAtomClass = get_basis_atom_class_from_ket(ket_1)
            basis = BasisAtomClass(ket_1.species, additional_kets=[ket_1, ket_2], database=self)
        else:
            if not system.is_diagonal:
                raise ValueError("System must be diagonal")
            basis = system.basis
        state_1, state_2 = (basis.get_corresponding_state(ket) for ket in (ket_1, ket_2))
        if system is not None:
            amplitudes = np.array(
                [state.get_amplitudes(ket)[0] for state, ket in zip((state_1, state_2), (ket_1, ket_2))]
            )
            if np.any(np.abs(amplitudes) ** 2 <= 0.5):
                raise ValueError("ket_1 or ket_2 does not clearly correspond to an eigenstate of the system.")
            if not np.allclose(np.abs(amplitudes), amplitudes):
                raise ValueError(
                    "The corresponding eigenstate of ket_1 or ket_2 has a non-positive amplitude, "
                    "this should not happen!"
                )

        cpp_operator_type = get_cpp_operator_type(operator)
        matrix_element_au = self._cpp.get_matrix_elements(state_1._cpp, state_2._cpp, cpp_operator_type, q)[0, 0]  # type: ignore
        matrix_element = QuantityScalar.from_base(matrix_element_au, operator)
        return matrix_element.to_unit(unit)


def get_basis_atom_class_from_ket(ket: "KetAtom") -> "type[BasisAtom]":
    import pairinteraction.backend._wrapped.basis.BasisAtom as BasisAtomModule

    type_ = type(ket._cpp).__name__.replace("KetAtom", "")
    try:
        return getattr(BasisAtomModule, f"BasisAtom{type_}")
    except AttributeError as err:
        raise ValueError(f"Unknown KetAtom {type(ket)}, cant find corresponding BasisAtom class") from err
