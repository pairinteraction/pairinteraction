from abc import ABC
from collections.abc import Sequence
from typing import TYPE_CHECKING, Union

import numpy as np

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.BasisCombined import (
    BasisCombinedBase,
    BasisCombinedComplexDouble,
    BasisCombinedComplexFloat,
    BasisCombinedDouble,
    BasisCombinedFloat,
)
from pairinteraction.unit_system import convert_base_to_quantity, convert_raw_to_base

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity, PlainUnit


class SystemCombinedBase(ABC):
    _CPPSystemCombined: Union[
        type[_backend.SystemCombinedFloat],
        type[_backend.SystemCombinedComplexFloat],
        type[_backend.SystemCombinedDouble],
        type[_backend.SystemCombinedComplexDouble],
    ]
    _BasisCombined: type[BasisCombinedBase]
    _DefaultDiagonalizer: Union[
        type[_backend.DiagonalizerEigenFloat],
        type[_backend.DiagonalizerEigenComplexFloat],
        type[_backend.DiagonalizerEigenDouble],
        type[_backend.DiagonalizerEigenComplexDouble],
    ]

    def __init__(self, basis: "BasisCombinedBase") -> None:
        self._cpp = self._CPPSystemCombined(basis._cpp)
        self.update_basis()

    @classmethod
    def _from_cpp_object(
        cls,
        cpp_obj: Union[
            _backend.SystemCombinedFloat,
            _backend.SystemCombinedComplexFloat,
            _backend.SystemCombinedDouble,
            _backend.SystemCombinedComplexDouble,
        ],
    ):
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        obj.update_basis()
        return obj

    @property
    def basis(self) -> "BasisCombinedBase":
        return self._basis

    def update_basis(self) -> None:
        self._basis = self._BasisCombined._from_cpp_object(self._cpp.get_basis())

    def set_order(self, order: int):
        self._cpp.set_order(order)
        return self

    def set_distance(self, distance: float, unit: Union[str, "PlainUnit"] = "mircrometer"):
        distance_base_units: float = convert_raw_to_base(distance, unit, "distance")
        self._cpp.set_distance(distance_base_units)
        return self

    def set_distance_vector(self, distance: Sequence[float], unit: Union[str, "PlainUnit"] = "mircrometer"):
        distance_base_units: np.ndarray = convert_raw_to_base(distance, unit, "distance")
        self._cpp.set_distance_vector(list(distance_base_units))
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

    def transformed(self, transformation) -> "SystemCombinedBase":
        new_cpp = self._cpp.transformed(transformation)
        return self._from_cpp_object(new_cpp)

    def get_eigenvalues(self, unit: Union[str, "PlainUnit"] = "GHz") -> np.ndarray:
        return self.get_eigenvalues_as_quantity(unit).magnitude

    def get_eigenvalues_as_quantity(self, unit: Union[str, "PlainUnit"] = "GHz") -> "PlainQuantity":
        eigenvalues_au = self._cpp.get_eigenvalues()
        return convert_base_to_quantity(eigenvalues_au, "energy", unit)

    def get_eigenbasis(self):
        return self._cpp.get_eigenbasis()


class SystemCombinedFloat(SystemCombinedBase):
    _cpp: _backend.SystemCombinedFloat
    _CPPSystemCombined = _backend.SystemCombinedFloat
    _BasisCombined = BasisCombinedFloat
    _DefaultDiagonalizer = _backend.DiagonalizerEigenFloat


class SystemCombinedComplexFloat(SystemCombinedBase):
    _cpp: _backend.SystemCombinedComplexFloat
    _CPPSystemCombined = _backend.SystemCombinedComplexFloat
    _BasisCombined = BasisCombinedComplexFloat
    _DefaultDiagonalizer = _backend.DiagonalizerEigenComplexFloat


class SystemCombinedDouble(SystemCombinedBase):
    _cpp: _backend.SystemCombinedDouble
    _CPPSystemCombined = _backend.SystemCombinedDouble
    _BasisCombined = BasisCombinedDouble
    _DefaultDiagonalizer = _backend.DiagonalizerEigenDouble


class SystemCombinedComplexDouble(SystemCombinedBase):
    _cpp: _backend.SystemCombinedComplexDouble
    _CPPSystemCombined = _backend.SystemCombinedComplexDouble
    _BasisCombined = BasisCombinedComplexDouble
    _DefaultDiagonalizer = _backend.DiagonalizerEigenComplexDouble
