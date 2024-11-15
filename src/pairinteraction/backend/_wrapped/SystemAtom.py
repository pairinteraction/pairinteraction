from abc import ABC
from collections.abc import Sequence
from typing import TYPE_CHECKING, Union

import numpy as np

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.BasisAtom import (
    BasisAtomBase,
    BasisAtomComplexDouble,
    BasisAtomComplexFloat,
    BasisAtomDouble,
    BasisAtomFloat,
)
from pairinteraction.unit_system import convert_base_to_quantity, convert_quantity_to_base, convert_raw_to_base

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity, PlainUnit


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

    def set_electric_field(self, electric_field: Sequence[float], unit: Union[str, "PlainUnit"] = "G"):
        electric_field_base_units: np.ndarray = convert_raw_to_base(electric_field, unit, "electric_field")
        self._cpp.set_electric_field(list(electric_field_base_units))
        return self

    def set_electric_field_from_quantity(self, electric_field: Sequence["PlainQuantity"]):
        electric_field_base_units: np.ndarray = convert_quantity_to_base(electric_field, "electric_field")
        self._cpp.set_electric_field(list(electric_field_base_units))
        return self

    def set_magnetic_field(self, magnetic_field: Sequence[float], unit: Union[str, "PlainUnit"] = "G"):
        magnetic_field_base_units: np.ndarray = convert_raw_to_base(magnetic_field, unit, "magnetic_field")
        self._cpp.set_magnetic_field(list(magnetic_field_base_units))
        return self

    def set_magnetic_field_from_quantity(self, magnetic_field: Sequence["PlainQuantity"]):
        magnetic_field_base_units: np.ndarray = convert_quantity_to_base(magnetic_field, "magnetic_field")
        self._cpp.set_magnetic_field(list(magnetic_field_base_units))
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

    def get_eigenvalues(self, unit: Union[str, "PlainUnit"] = "GHz") -> np.ndarray:
        return self.get_eigenvalues_as_quantity(unit).magnitude

    def get_eigenvalues_as_quantity(self, unit: Union[str, "PlainUnit"] = "GHz") -> "PlainQuantity":
        eigenvalues_au = self._cpp.get_eigenvalues()
        return convert_base_to_quantity(eigenvalues_au, "energy", unit)

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
