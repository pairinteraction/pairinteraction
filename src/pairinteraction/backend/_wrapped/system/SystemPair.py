from collections.abc import Collection
from typing import TYPE_CHECKING, Any, ClassVar, TypeVar, Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.basis.BasisPair import (
    BasisPairComplexDouble,
    BasisPairComplexFloat,
    BasisPairDouble,
    BasisPairFloat,
)
from pairinteraction.backend._wrapped.system.System import SystemBase
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

    from pairinteraction.units import Array

    SelfSystemPair_t = TypeVar("SelfSystemPair_t", bound="SystemPair")

Basis_t = TypeVar(
    "Basis_t", "BasisPairFloat", "BasisPairComplexFloat", "BasisPairDouble", "BasisPairComplexDouble"
)
UnionCPPSystemPair = Union[
    _backend.SystemPairFloat,
    _backend.SystemPairComplexFloat,
    _backend.SystemPairDouble,
    _backend.SystemPairComplexDouble,
]
UnionTypeCPPSystemPair = Union[
    type[_backend.SystemPairFloat],
    type[_backend.SystemPairComplexFloat],
    type[_backend.SystemPairDouble],
    type[_backend.SystemPairComplexDouble],
]


class SystemPairBase(SystemBase[Basis_t]):
    _cpp: UnionCPPSystemPair
    _cpp_type: ClassVar[UnionTypeCPPSystemPair]

    def set_order(self: "SelfSystemPair_t", order: int) -> "SelfSystemPair_t":
        self._cpp.set_order(order)
        return self

    def set_distance(
        self: "SelfSystemPair_t", distance: Union[float, "PlainQuantity[float]"], unit: str = "pint"
    ) -> "SelfSystemPair_t":
        distance_au = QuantityScalar(distance, unit).to_base("DISTANCE")
        self._cpp.set_distance(distance_au)
        return self

    def set_distance_vector(
        self: "SelfSystemPair_t",
        distance: Union["PlainQuantity[Array]", Collection[Union[float, "PlainQuantity[float]"]]],
        unit: str = "pint",
    ) -> "SelfSystemPair_t":
        distance_au = [QuantityScalar(v, unit).to_base("DISTANCE") for v in distance]
        self._cpp.set_distance_vector(distance_au)
        return self


class SystemPairFloat(SystemPairBase[BasisPairFloat]):
    _cpp: _backend.SystemPairFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemPairFloat
    _TypeBasis = BasisPairFloat


class SystemPairComplexFloat(SystemPairBase[BasisPairComplexFloat]):
    _cpp: _backend.SystemPairComplexFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemPairComplexFloat
    _TypeBasis = BasisPairComplexFloat


class SystemPairDouble(SystemPairBase[BasisPairDouble]):
    _cpp: _backend.SystemPairDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemPairDouble
    _TypeBasis = BasisPairDouble


class SystemPairComplexDouble(SystemPairBase[BasisPairComplexDouble]):
    _cpp: _backend.SystemPairComplexDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemPairComplexDouble
    _TypeBasis = BasisPairComplexDouble


SystemPair = SystemPairBase[Any]
