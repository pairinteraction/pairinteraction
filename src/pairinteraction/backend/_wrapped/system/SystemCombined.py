from collections.abc import Collection
from typing import TYPE_CHECKING, Any, ClassVar, TypeVar, Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.basis.BasisCombined import (
    BasisCombinedComplexDouble,
    BasisCombinedComplexFloat,
    BasisCombinedDouble,
    BasisCombinedFloat,
)
from pairinteraction.backend._wrapped.system.System import SystemBase
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

    from pairinteraction.units import Array

    SelfSystemCombined_t = TypeVar("SelfSystemCombined_t", bound="SystemCombined")

Basis_t = TypeVar(
    "Basis_t", "BasisCombinedFloat", "BasisCombinedComplexFloat", "BasisCombinedDouble", "BasisCombinedComplexDouble"
)
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


class SystemCombinedBase(SystemBase[Basis_t]):
    _cpp: UnionCPPSystemCombined
    _cpp_type: ClassVar[UnionTypeCPPSystemCombined]

    def set_order(self: "SelfSystemCombined_t", order: int) -> "SelfSystemCombined_t":
        self._cpp.set_order(order)
        return self

    def set_distance(
        self: "SelfSystemCombined_t", distance: Union[float, "PlainQuantity[float]"], unit: str = "pint"
    ) -> "SelfSystemCombined_t":
        distance_au = QuantityScalar(distance, unit).to_base("DISTANCE")
        self._cpp.set_distance(distance_au)
        return self

    def set_distance_vector(
        self: "SelfSystemCombined_t",
        distance: Union["PlainQuantity[Array]", Collection[Union[float, "PlainQuantity[float]"]]],
        unit: str = "pint",
    ) -> "SelfSystemCombined_t":
        distance_au = [QuantityScalar(v, unit).to_base("DISTANCE") for v in distance]
        self._cpp.set_distance_vector(distance_au)
        return self


class SystemCombinedFloat(SystemCombinedBase[BasisCombinedFloat]):
    _cpp: _backend.SystemCombinedFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemCombinedFloat
    _TypeBasis = BasisCombinedFloat


class SystemCombinedComplexFloat(SystemCombinedBase[BasisCombinedComplexFloat]):
    _cpp: _backend.SystemCombinedComplexFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemCombinedComplexFloat
    _TypeBasis = BasisCombinedComplexFloat


class SystemCombinedDouble(SystemCombinedBase[BasisCombinedDouble]):
    _cpp: _backend.SystemCombinedDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemCombinedDouble
    _TypeBasis = BasisCombinedDouble


class SystemCombinedComplexDouble(SystemCombinedBase[BasisCombinedComplexDouble]):
    _cpp: _backend.SystemCombinedComplexDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemCombinedComplexDouble
    _TypeBasis = BasisCombinedComplexDouble


SystemCombined = SystemCombinedBase[Any]
