from abc import ABC
from typing import TYPE_CHECKING, Any, ClassVar, Generic, Optional, TypeVar, Union, overload

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.Diagonalizer import Diagonalizer
from pairinteraction.backend._wrapped.get_functions import get_cpp_diagonalizer
from pairinteraction.units import QuantityArray, QuantityScalar, QuantitySparse

if TYPE_CHECKING:
    import numpy as np
    import numpy.typing as npt
    from pint.facets.plain import PlainQuantity
    from scipy.sparse import csr_matrix
    from typing_extensions import Self

    from pairinteraction.backend._wrapped.basis.BasisAtom import BasisAtom
    from pairinteraction.backend._wrapped.basis.BasisPair import BasisPair
    from pairinteraction.units import Array


BasisType = TypeVar("BasisType", bound=Union["BasisAtom", "BasisPair"])
UnionCPPSystem = Any
# UnionCPPSystem is supposed to be System(|System)(Atom|Pair)(Real|Complex)
UnionTypeCPPSystem = Any


class SystemBase(ABC, Generic[BasisType]):
    """Base class for all System objects.

    The system objects are meant to represent the full physical system, including the basis
    as well as possible fields and interactions which form the Hamiltonian.

    Typically, the system objects are simply constructed from a basis object.
    The fields and interactions can then be set afterwards.

    All system objects share a few common attributes and methods, that are defined in this base class, e.g.:
    - the basis of the system,
    - the matrix stored as scipy sparse matrix, which represents the Hamiltonian,
    - ...
    """

    _cpp: UnionCPPSystem
    _cpp_type: ClassVar[UnionTypeCPPSystem]
    _basis: BasisType
    _TypeBasis: type[BasisType]  # should by ClassVar, but cannot be nested yet

    def __init__(self, basis: BasisType) -> None:
        self._cpp = self._cpp_type(basis._cpp)  # type: ignore [reportPrivateUsage]
        self.update_basis()

    @classmethod
    def _from_cpp_object(cls: "type[Self]", cpp_obj: UnionCPPSystem) -> "Self":
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        obj.update_basis()
        return obj

    def __str__(self) -> str:
        return self.__repr__()

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.basis!r}, is_diagonal={self.is_diagonal})"

    def update_basis(self) -> None:
        self._basis = self._TypeBasis._from_cpp_object(self._cpp.get_basis())  # type: ignore

    def diagonalize(
        self,
        diagonalizer: Diagonalizer = "Eigen",
        precision: int = 12,
        sort_by_energy: bool = True,
        energy_range: Union[tuple[float, float], tuple["PlainQuantity[float]", "PlainQuantity[float]"], None] = None,
        energy_unit: Optional[str] = None,
    ) -> "Self":
        cpp_diagonalizer = get_cpp_diagonalizer(diagonalizer, self._cpp)
        if energy_range is None:
            self._cpp.diagonalize(cpp_diagonalizer, precision)
        else:
            min_energy_au = QuantityScalar(energy_range[0], energy_unit).to_base("ENERGY")
            max_energy_au = QuantityScalar(energy_range[1], energy_unit).to_base("ENERGY")
            cpp_range = _backend.RangeDouble(min_energy_au, max_energy_au)
            self._cpp.diagonalize(cpp_diagonalizer, precision, cpp_range)
        if sort_by_energy:
            sorter = self._cpp.get_sorter([_backend.TransformationType.SORT_BY_ENERGY])
            self._cpp.transform(sorter)

        self.update_basis()
        return self

    @property
    def is_diagonal(self) -> bool:
        return self._cpp.is_diagonal()

    @property
    def basis(self) -> BasisType:
        return self._basis

    @property
    def matrix(self) -> "csr_matrix":
        return self._cpp.get_matrix()

    def get_eigenbasis(self) -> BasisType:
        cpp_eigenbasis = self._cpp.get_eigenbasis()
        return self._TypeBasis._from_cpp_object(cpp_eigenbasis)  # type: ignore

    @overload
    def get_eigenvalues(self) -> "PlainQuantity[Array]": ...

    @overload
    def get_eigenvalues(self, unit: str) -> "npt.NDArray[np.floating[Any]]": ...

    def get_eigenvalues(self, unit: Optional[str] = None):
        eigenvalues_au = self._cpp.get_eigenvalues()
        eigenvalues = QuantityArray.from_base(eigenvalues_au, "ENERGY")
        return eigenvalues.to_unit(unit)

    @overload
    def get_hamiltonian(self) -> "PlainQuantity[csr_matrix]": ...  # type: ignore [reportInvalidTypeArguments]

    @overload
    def get_hamiltonian(self, unit: str) -> "csr_matrix": ...

    def get_hamiltonian(self, unit: Optional[str] = None):
        hamiltonian_au = self._cpp.get_matrix()
        hamiltonian = QuantitySparse.from_base(hamiltonian_au, "ENERGY")
        return hamiltonian.to_unit(unit)


System = SystemBase[Any]
