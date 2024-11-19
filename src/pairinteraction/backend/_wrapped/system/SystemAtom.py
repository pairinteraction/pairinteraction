from collections.abc import Collection
from typing import TYPE_CHECKING, Union

import pairinteraction.backend._backend as _backend
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
    _cpp_type: UnionTypeCPPSystemAtom
    _TypeBasis: type[BasisAtomBase]

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


class SystemAtomFloat(SystemAtomBase):
    _cpp: _backend.SystemAtomFloat
    _cpp_type = _backend.SystemAtomFloat
    _TypeBasis = BasisAtomFloat


class SystemAtomComplexFloat(SystemAtomBase):
    _cpp: _backend.SystemAtomComplexFloat
    _cpp_type = _backend.SystemAtomComplexFloat
    _TypeBasis = BasisAtomComplexFloat


class SystemAtomDouble(SystemAtomBase):
    _cpp: _backend.SystemAtomDouble
    _cpp_type = _backend.SystemAtomDouble
    _TypeBasis = BasisAtomDouble


class SystemAtomComplexDouble(SystemAtomBase):
    _cpp: _backend.SystemAtomComplexDouble
    _cpp_type = _backend.SystemAtomComplexDouble
    _TypeBasis = BasisAtomComplexDouble
