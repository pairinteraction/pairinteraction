from abc import ABC
from typing import TYPE_CHECKING, Any, Generic, Optional, TypeVar, Union, overload

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.basis.Basis import BasisBase
from pairinteraction.backend._wrapped.Diagonalizer import Diagonalizer, get_cpp_diagonalizer
from pairinteraction.unit_system import Qties

if TYPE_CHECKING:
    from numpy.typing import ArrayLike
    from pint.facets.plain import PlainQuantity

    Self = TypeVar("Self", bound="SystemBase")

# UnionCPPSystem = Union[
#     _backend.SystemSystem...,
# ]
UnionCPPSystem = Any

# UnionTypeCPPSystem = Union[
#     type[_backend.SystemAtomCreatorFloat],
# ]
UnionTypeCPPSystem = Any

UnionCPPRange = Union[_backend.RangeFloat, _backend.RangeDouble]

Basis_t = TypeVar("Basis_t", bound="BasisBase")


class SystemBase(ABC, Generic[Basis_t]):
    _cpp: UnionCPPSystem
    _cpp_type: UnionTypeCPPSystem
    _TypeBasis: type[BasisBase]

    def __init__(self, basis: Basis_t) -> None:
        self._cpp = self._cpp_type(basis._cpp)
        self.update_basis()

    @classmethod
    def _from_cpp_object(cls, cpp_obj: UnionCPPSystem):
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        obj.update_basis()
        return obj

    def update_basis(self) -> None:
        self._basis: Basis_t = self._TypeBasis._from_cpp_object(self._cpp.get_basis())

    def diagonalize(
        self,
        diagonalizer: Diagonalizer = "Eigen",
        precision: int = 12,
        eigenvalue_range: Optional[UnionCPPRange] = None,
    ):
        cpp_diagonalizer = get_cpp_diagonalizer(diagonalizer, self._cpp)
        if eigenvalue_range is None:
            self._cpp.diagonalize(cpp_diagonalizer, precision)
        else:
            self._cpp.diagonalize(cpp_diagonalizer, precision, eigenvalue_range)
        self.update_basis()
        return self

    @property
    def basis(self) -> Basis_t:
        return self._basis

    @property
    def matrix(self):
        return self._cpp.get_matrix()

    @overload
    def get_eigenvalues(self) -> "PlainQuantity[float]": ...

    @overload
    def get_eigenvalues(self, unit: str) -> "ArrayLike": ...

    def get_eigenvalues(self, unit: str = "pint"):
        eigenvalues_au = self._cpp.get_eigenvalues()
        eigenvalues = Qties.from_base(eigenvalues_au, "energy")
        return eigenvalues.to_unit(unit)

    def get_eigenbasis(self) -> Basis_t:
        cpp_eigenbasis = self._cpp.get_eigenbasis()
        return self._TypeBasis._from_cpp_object(cpp_eigenbasis)

    def transform(self, transformation) -> "Self":
        self._cpp.transform(transformation)
        return self
