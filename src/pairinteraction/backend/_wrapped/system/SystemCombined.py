from collections.abc import Collection
from typing import TYPE_CHECKING, Union

import pairinteraction.backend._backend as _backend
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


class SystemCombinedBase(SystemBase):
    _cpp: UnionCPPSystemCombined
    _cpp_type: UnionTypeCPPSystemCombined
    _TypeBasis: type[BasisCombinedBase]

    def set_order(self, order: int):
        self._cpp.set_order(order)
        return self

    def set_distance(self, distance: Union[float, "PlainQuantity"], unit: str = "pint"):
        distance_au = Qty(distance, unit).to_base("distance")
        self._cpp.set_distance(distance_au)
        return self

    def set_distance_vector(
        self, distance: Union["PlainQuantity", Collection[Union[float, "PlainQuantity"]]], unit: str = "pint"
    ):
        distance_au = [Qty(v, unit).to_base("distance") for v in distance]
        self._cpp.set_distance_vector(distance_au)
        return self


class SystemCombinedFloat(SystemCombinedBase):
    _cpp: _backend.SystemCombinedFloat
    _cpp_type = _backend.SystemCombinedFloat
    _TypeBasis = BasisCombinedFloat


class SystemCombinedComplexFloat(SystemCombinedBase):
    _cpp: _backend.SystemCombinedComplexFloat
    _cpp_type = _backend.SystemCombinedComplexFloat
    _TypeBasis = BasisCombinedComplexFloat


class SystemCombinedDouble(SystemCombinedBase):
    _cpp: _backend.SystemCombinedDouble
    _cpp_type = _backend.SystemCombinedDouble
    _TypeBasis = BasisCombinedDouble


class SystemCombinedComplexDouble(SystemCombinedBase):
    _cpp: _backend.SystemCombinedComplexDouble
    _cpp_type = _backend.SystemCombinedComplexDouble
    _TypeBasis = BasisCombinedComplexDouble
