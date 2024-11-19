from typing import TYPE_CHECKING, Optional, Union

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._backend import (
    KetCombinedComplexDouble,
    KetCombinedComplexFloat,
    KetCombinedDouble,
    KetCombinedFloat,
)
from pairinteraction.backend._wrapped.Basis import BasisBase
from pairinteraction.backend._wrapped.SystemAtom import SystemAtomBase
from pairinteraction.unit_system import Qty

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

UnionCPPBasisCombined = Union[
    _backend.BasisCombinedFloat,
    _backend.BasisCombinedComplexFloat,
    _backend.BasisCombinedDouble,
    _backend.BasisCombinedComplexDouble,
]
UnionTypeCPPBasisCombinedCreator = Union[
    type[_backend.BasisCombinedCreatorFloat],
    type[_backend.BasisCombinedCreatorComplexFloat],
    type[_backend.BasisCombinedCreatorDouble],
    type[_backend.BasisCombinedCreatorComplexDouble],
]
UnionTypeKetCombined = Union[
    type[KetCombinedFloat],
    type[KetCombinedComplexFloat],
    type[KetCombinedDouble],
    type[KetCombinedComplexDouble],
]
UnionKetCombined = Union[
    KetCombinedFloat,
    KetCombinedComplexFloat,
    KetCombinedDouble,
    KetCombinedComplexDouble,
]


class BasisCombinedBase(BasisBase):
    _cpp: UnionCPPBasisCombined
    _cpp_creator: UnionTypeCPPBasisCombinedCreator
    _TypeKet: UnionTypeKetCombined

    def __init__(
        self,
        systems: list[SystemAtomBase],
        m: Optional[tuple[float, float]] = None,
        energy: Union[tuple[float, float], tuple["PlainQuantity", "PlainQuantity"], None] = None,
        energy_unit: str = "pint",
    ) -> None:
        creator = self._cpp_creator()
        for system in systems:
            creator.add(system._cpp)
        if m is not None:
            creator.restrict_quantum_number_m(*m)
        if energy is not None:
            min_energy_au = Qty(energy[0], energy_unit).to_base("energy")
            max_energy_au = Qty(energy[1], energy_unit).to_base("energy")
            creator.restrict_energy(min_energy_au, max_energy_au)
        self._cpp = creator.create()


class BasisCombinedFloat(BasisCombinedBase):
    _cpp: _backend.BasisCombinedFloat
    _cpp_creator = _backend.BasisCombinedCreatorFloat
    _TypeKet = KetCombinedFloat


class BasisCombinedComplexFloat(BasisCombinedBase):
    _cpp: _backend.BasisCombinedComplexFloat
    _cpp_creator = _backend.BasisCombinedCreatorComplexFloat
    _TypeKet = KetCombinedComplexFloat


class BasisCombinedDouble(BasisCombinedBase):
    _cpp: _backend.BasisCombinedDouble
    _cpp_creator = _backend.BasisCombinedCreatorDouble
    _TypeKet = KetCombinedDouble


class BasisCombinedComplexDouble(BasisCombinedBase):
    _cpp: _backend.BasisCombinedComplexDouble
    _cpp_creator = _backend.BasisCombinedCreatorComplexDouble
    _TypeKet = KetCombinedComplexDouble
