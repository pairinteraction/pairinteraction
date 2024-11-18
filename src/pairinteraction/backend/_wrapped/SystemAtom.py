from abc import ABC
from collections.abc import Collection
from typing import TYPE_CHECKING, Union

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.BasisAtom import (
    BasisAtomBase,
    BasisAtomComplexDouble,
    BasisAtomComplexFloat,
    BasisAtomDouble,
    BasisAtomFloat,
)
from pairinteraction.unit_system import Qties, Qty

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity


class SystemAtomBase(ABC):
    _CPPSystemAtom: Union[
        type[_backend.SystemAtomFloat],
        type[_backend.SystemAtomComplexFloat],
        type[_backend.SystemAtomDouble],
        type[_backend.SystemAtomComplexDouble],
    ]
    _BasisAtom: type[BasisAtomBase]
    _DefaultDiagonalizer: Union[
        type[_backend.DiagonalizerEigenFloat],
        type[_backend.DiagonalizerEigenComplexFloat],
        type[_backend.DiagonalizerEigenDouble],
        type[_backend.DiagonalizerEigenComplexDouble],
    ]

    def __init__(self, basis: "BasisAtomBase") -> None:
        self._cpp = self._CPPSystemAtom(basis._cpp)
        self.update_basis()

    @classmethod
    def _from_cpp_object(
        cls,
        cpp_obj: Union[
            _backend.SystemAtomFloat,
            _backend.SystemAtomComplexFloat,
            _backend.SystemAtomDouble,
            _backend.SystemAtomComplexDouble,
        ],
    ):
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        obj.update_basis()
        return obj

    @property
    def basis(self) -> "BasisAtomBase":
        return self._basis

    def update_basis(self) -> None:
        self._basis = self._BasisAtom._from_cpp_object(self._cpp.get_basis())

    def set_electric_field(
        self, electric_field: Union["PlainQuantity", Collection[Union[float, "PlainQuantity"]]], unit: str = "pint"
    ):
        electric_field_au = [Qty(v, unit).to_base("electric_field") for v in electric_field]
        self._cpp.set_electric_field(electric_field_au)
        return self

    def set_magnetic_field(
        self, magnetic_field: Union["PlainQuantity", Collection[Union[float, "PlainQuantity"]]], unit: str = "pint"
    ):
        magnetic_field_au = [Qty(v, unit).to_base("magnetic_field") for v in magnetic_field]
        self._cpp.set_magnetic_field(magnetic_field_au)
        return self

    def enable_diamagnetism(self, enable: bool):
        self._cpp.enable_diamagnetism(enable)
        return self

    def diagonalize(
        self,
        diagonalizer: Union[_backend.DiagonalizerInterfaceFloat, _backend.DiagonalizerInterfaceDouble, None] = None,
        precision: int = 12,
        eigenvalue_range: Union[_backend.RangeFloat, _backend.RangeDouble, None] = None,
    ):
        if diagonalizer is None:
            diagonalizer = self._DefaultDiagonalizer()
        if eigenvalue_range is None:
            self._cpp.diagonalize(diagonalizer, precision)
        else:
            self._cpp.diagonalize(diagonalizer, precision, eigenvalue_range)
        self.update_basis()
        return self

    @property
    def matrix(self):
        return self._cpp.get_matrix()

    def transform(self, transformation) -> "SystemAtomBase":
        self._cpp.transform(transformation)
        return self

    def get_eigenvalues(self, unit: str = "pint"):
        eigenvalues_au = self._cpp.get_eigenvalues()
        eigenvalues = Qties.from_base(eigenvalues_au, "energy")
        return eigenvalues.to_unit(unit)

    def get_eigenbasis(self):
        return self._cpp.get_eigenbasis()


class SystemAtomFloat(SystemAtomBase):
    _cpp: _backend.SystemAtomFloat
    _CPPSystemAtom = _backend.SystemAtomFloat
    _BasisAtom = BasisAtomFloat
    _DefaultDiagonalizer = _backend.DiagonalizerEigenFloat


class SystemAtomComplexFloat(SystemAtomBase):
    _cpp: _backend.SystemAtomComplexFloat
    _CPPSystemAtom = _backend.SystemAtomComplexFloat
    _BasisAtom = BasisAtomComplexFloat
    _DefaultDiagonalizer = _backend.DiagonalizerEigenComplexFloat


class SystemAtomDouble(SystemAtomBase):
    _cpp: _backend.SystemAtomDouble
    _CPPSystemAtom = _backend.SystemAtomDouble
    _BasisAtom = BasisAtomDouble
    _DefaultDiagonalizer = _backend.DiagonalizerEigenDouble


class SystemAtomComplexDouble(SystemAtomBase):
    _cpp: _backend.SystemAtomComplexDouble
    _CPPSystemAtom = _backend.SystemAtomComplexDouble
    _BasisAtom = BasisAtomComplexDouble
    _DefaultDiagonalizer = _backend.DiagonalizerEigenComplexDouble
