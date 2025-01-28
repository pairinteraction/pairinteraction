from collections.abc import Collection
from typing import TYPE_CHECKING, Any, ClassVar, TypeVar, Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.basis.BasisAtom import (
    BasisAtomComplexDouble,
    BasisAtomComplexFloat,
    BasisAtomDouble,
    BasisAtomFloat,
)
from pairinteraction.backend._wrapped.system.System import SystemBase
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity
    from typing_extensions import Self

    from pairinteraction.units import Array

Basis_t = TypeVar("Basis_t", "BasisAtomFloat", "BasisAtomComplexFloat", "BasisAtomDouble", "BasisAtomComplexDouble")
UnionCPPSystemAtom = Union[
    _backend.SystemAtomFloat,
    _backend.SystemAtomComplexFloat,
    _backend.SystemAtomDouble,
    _backend.SystemAtomComplexDouble,
]
UnionTypeCPPSystemAtom = Union[
    type[_backend.SystemAtomFloat],
    type[_backend.SystemAtomComplexFloat],
    type[_backend.SystemAtomDouble],
    type[_backend.SystemAtomComplexDouble],
]


class SystemAtomBase(SystemBase[Basis_t]):
    _cpp: UnionCPPSystemAtom
    _cpp_type: ClassVar[UnionTypeCPPSystemAtom]

    def set_electric_field(
        self: "Self",
        electric_field: Union["PlainQuantity[Array]", Collection[Union[float, "PlainQuantity[float]"]]],
        unit: str = "pint",
    ) -> "Self":
        electric_field_au = [QuantityScalar(v, unit).to_base("ELECTRIC_FIELD") for v in electric_field]
        self._cpp.set_electric_field(electric_field_au)
        return self

    def set_magnetic_field(
        self: "Self",
        magnetic_field: Union["PlainQuantity[Array]", Collection[Union[float, "PlainQuantity[float]"]]],
        unit: str = "pint",
    ) -> "Self":
        magnetic_field_au = [QuantityScalar(v, unit).to_base("MAGNETIC_FIELD") for v in magnetic_field]
        self._cpp.set_magnetic_field(magnetic_field_au)
        return self

    def enable_diamagnetism(self: "Self", enable: bool = True) -> "Self":
        self._cpp.enable_diamagnetism(enable)
        return self


class SystemAtomFloat(SystemAtomBase[BasisAtomFloat]):
    _cpp: _backend.SystemAtomFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemAtomFloat
    _TypeBasis = BasisAtomFloat


class SystemAtomComplexFloat(SystemAtomBase[BasisAtomComplexFloat]):
    _cpp: _backend.SystemAtomComplexFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemAtomComplexFloat
    _TypeBasis = BasisAtomComplexFloat


class SystemAtomDouble(SystemAtomBase[BasisAtomDouble]):
    _cpp: _backend.SystemAtomDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemAtomDouble
    _TypeBasis = BasisAtomDouble


class SystemAtomComplexDouble(SystemAtomBase[BasisAtomComplexDouble]):
    _cpp: _backend.SystemAtomComplexDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemAtomComplexDouble
    _TypeBasis = BasisAtomComplexDouble


SystemAtom = SystemAtomBase[Any]
