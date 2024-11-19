from collections.abc import Collection
from typing import TYPE_CHECKING, ClassVar, Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.basis.BasisAtom import (
    BasisAtomBase,
    BasisAtomComplexDouble,
    BasisAtomComplexFloat,
    BasisAtomDouble,
    BasisAtomFloat,
)
from pairinteraction.backend._wrapped.system.System import SystemBase
from pairinteraction.unit_system import Qty

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

    from pairinteraction.unit_system import Array

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


class SystemAtomBase(SystemBase[BasisAtomBase]):
    _cpp: UnionCPPSystemAtom
    _cpp_type: ClassVar[UnionTypeCPPSystemAtom]

    def set_electric_field(
        self,
        electric_field: Union["PlainQuantity[Array]", Collection[Union[float, "PlainQuantity[float]"]]],
        unit: str = "pint",
    ):
        electric_field_au = [Qty(v, unit).to_base("electric_field") for v in electric_field]
        self._cpp.set_electric_field(electric_field_au)
        return self

    def set_magnetic_field(
        self,
        magnetic_field: Union["PlainQuantity[Array]", Collection[Union[float, "PlainQuantity[float]"]]],
        unit: str = "pint",
    ):
        magnetic_field_au = [Qty(v, unit).to_base("magnetic_field") for v in magnetic_field]
        self._cpp.set_magnetic_field(magnetic_field_au)
        return self

    def enable_diamagnetism(self, enable: bool = True):
        self._cpp.enable_diamagnetism(enable)
        return self


class SystemAtomFloat(SystemAtomBase):
    _cpp: _backend.SystemAtomFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemAtomFloat
    _TypeBasis = BasisAtomFloat


class SystemAtomComplexFloat(SystemAtomBase):
    _cpp: _backend.SystemAtomComplexFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemAtomComplexFloat
    _TypeBasis = BasisAtomComplexFloat


class SystemAtomDouble(SystemAtomBase):
    _cpp: _backend.SystemAtomDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemAtomDouble
    _TypeBasis = BasisAtomDouble


class SystemAtomComplexDouble(SystemAtomBase):
    _cpp: _backend.SystemAtomComplexDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemAtomComplexDouble
    _TypeBasis = BasisAtomComplexDouble
