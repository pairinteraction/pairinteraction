from abc import ABC
from typing import TYPE_CHECKING, Any, ClassVar, Generic, TypeVar, Union, overload

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.Diagonalizer import Diagonalizer, get_cpp_diagonalizer, get_cpp_range
from pairinteraction.units import QuantityArray, QuantityScalar, QuantitySparse

if TYPE_CHECKING:
    import numpy as np
    import numpy.typing as npt
    from pint.facets.plain import PlainQuantity
    from scipy.sparse import csr_matrix

    from pairinteraction.backend._wrapped.basis.BasisAtom import BasisAtom
    from pairinteraction.backend._wrapped.basis.BasisPair import BasisPair
    from pairinteraction.units import Array

    SelfSystem_t = TypeVar("SelfSystem_t", bound="System")

Basis_t = TypeVar("Basis_t", bound=Union["BasisAtom", "BasisPair"])
UnionCPPSystem = Any
# UnionCPPSystem is supposed to be System(|System)(Atom|Pair)(Float|Double|ComplexFloat|ComplexDouble)
UnionTypeCPPSystem = Any
# UnionTypeCPPSystem is supposed to be type[System(Atom|Pair)(Float|Double|ComplexFloat|ComplexDouble)]
UnionCPPRange = Union[_backend.RangeFloat, _backend.RangeDouble]


class SystemBase(ABC, Generic[Basis_t]):
    _cpp: UnionCPPSystem
    _cpp_type: ClassVar[UnionTypeCPPSystem]
    _basis: Basis_t
    _TypeBasis: type[Basis_t]  # should by ClassVar, but cannot be nested yet

    def __init__(self, basis: Basis_t) -> None:
        self._cpp = self._cpp_type(basis._cpp)  # type: ignore [reportPrivateUsage]
        self.update_basis()

    @classmethod
    def _from_cpp_object(cls: "type[SelfSystem_t]", cpp_obj: UnionCPPSystem) -> "SelfSystem_t":
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        obj.update_basis()
        return obj

    def __str__(self) -> str:
        return self.__repr__()

    def __repr__(self) -> str:
        return f"{type(self).__name__}({repr(self.basis)}, is_diagonal={self.is_diagonal})"

    def update_basis(self) -> None:
        self._basis = self._TypeBasis._from_cpp_object(self._cpp.get_basis())  # type: ignore

    def diagonalize(
        self,
        diagonalizer: Diagonalizer = "Eigen",
        precision: int = 12,
        sort_by_energy: bool = True,
        energy_range: Union[tuple[float, float], tuple["PlainQuantity[float]", "PlainQuantity[float]"], None] = None,
        energy_unit: str = "pint",
    ):
        cpp_diagonalizer = get_cpp_diagonalizer(diagonalizer, self._cpp)
        if energy_range is None:
            self._cpp.diagonalize(cpp_diagonalizer, precision)
        else:
            min_energy_au = QuantityScalar(energy_range[0], energy_unit).to_base("ENERGY")
            max_energy_au = QuantityScalar(energy_range[1], energy_unit).to_base("ENERGY")
            cpp_range_class = get_cpp_range(self._cpp)
            cpp_range = cpp_range_class(min_energy_au, max_energy_au)
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
    def basis(self) -> Basis_t:
        return self._basis

    @property
    def matrix(self):
        return self._cpp.get_matrix()

    @overload
    def get_eigenvalues(self) -> "PlainQuantity[Array]": ...

    @overload
    def get_eigenvalues(self, unit: str) -> "npt.NDArray[np.floating[Any]]": ...

    def get_eigenvalues(self, unit: str = "pint"):
        eigenvalues_au = self._cpp.get_eigenvalues()
        eigenvalues = QuantityArray.from_base(eigenvalues_au, "ENERGY")
        return eigenvalues.to_unit(unit)

    def get_eigenbasis(self) -> Basis_t:
        cpp_eigenbasis = self._cpp.get_eigenbasis()
        return self._TypeBasis._from_cpp_object(cpp_eigenbasis)  # type: ignore

    @overload
    def get_hamiltonian(self) -> "PlainQuantity[csr_matrix]": ...  # type: ignore [reportInvalidTypeArguments]

    @overload
    def get_hamiltonian(self, unit: str) -> "csr_matrix": ...

    def get_hamiltonian(self, unit: str = "pint"):
        hamiltonian_au = self._cpp.get_matrix()
        hamiltonian = QuantitySparse.from_base(hamiltonian_au, "ENERGY")
        return hamiltonian.to_unit(unit)


System = SystemBase[Any]
