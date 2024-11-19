from collections.abc import Collection
from typing import TYPE_CHECKING, ClassVar, Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.basis.BasisCombined import (
    BasisCombinedBase,
    BasisCombinedComplexDouble,
    BasisCombinedComplexFloat,
    BasisCombinedDouble,
    BasisCombinedFloat,
)
from pairinteraction.backend._wrapped.system.System import SystemBase
from pairinteraction.unit_system import Qty

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

    from pairinteraction.unit_system import Array

UnionCPPSystemCombined = Union[
    _backend.SystemCombinedFloat,
    _backend.SystemCombinedComplexFloat,
    _backend.SystemCombinedDouble,
    _backend.SystemCombinedComplexDouble,
]
UnionTypeCPPSystemCombined = Union[
    type[_backend.SystemCombinedFloat],
    type[_backend.SystemCombinedComplexFloat],
    type[_backend.SystemCombinedDouble],
    type[_backend.SystemCombinedComplexDouble],
]


class SystemCombinedBase(SystemBase[BasisCombinedBase]):
    _cpp: UnionCPPSystemCombined
    _cpp_type: ClassVar[UnionTypeCPPSystemCombined]

    def set_order(self, order: int):
        self._cpp.set_order(order)
        return self

    def set_distance(self, distance: Union[float, "PlainQuantity[float]"], unit: str = "pint"):
        distance_au = Qty(distance, unit).to_base("distance")
        self._cpp.set_distance(distance_au)
        return self

    def set_distance_vector(
        self,
        distance: Union["PlainQuantity[Array]", Collection[Union[float, "PlainQuantity[float]"]]],
        unit: str = "pint",
    ):
        distance_au = [Qty(v, unit).to_base("distance") for v in distance]
        self._cpp.set_distance_vector(distance_au)
        return self


class SystemCombinedFloat(SystemCombinedBase):
    _cpp: _backend.SystemCombinedFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemCombinedFloat
    _TypeBasis = BasisCombinedFloat


class SystemCombinedComplexFloat(SystemCombinedBase):
    _cpp: _backend.SystemCombinedComplexFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemCombinedComplexFloat
    _TypeBasis = BasisCombinedComplexFloat


class SystemCombinedDouble(SystemCombinedBase):
    _cpp: _backend.SystemCombinedDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemCombinedDouble
    _TypeBasis = BasisCombinedDouble


class SystemCombinedComplexDouble(SystemCombinedBase):
    _cpp: _backend.SystemCombinedComplexDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemCombinedComplexDouble
    _TypeBasis = BasisCombinedComplexDouble
