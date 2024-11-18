from functools import cached_property
from typing import TYPE_CHECKING, Optional, Union

import scipy.sparse

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._backend import (
    KetCombinedComplexDouble,
    KetCombinedComplexFloat,
    KetCombinedDouble,
    KetCombinedFloat,
)
from pairinteraction.backend._wrapped.SystemAtom import SystemAtomBase
from pairinteraction.unit_system import Qty

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity


class BasisCombinedBase:
    _BasisCombinedCreator: Union[
        type[_backend.BasisCombinedCreatorFloat],
        type[_backend.BasisCombinedCreatorComplexFloat],
        type[_backend.BasisCombinedCreatorDouble],
        type[_backend.BasisCombinedCreatorComplexDouble],
    ]
    _KetCombined: Union[
        type[KetCombinedFloat], type[KetCombinedComplexFloat], type[KetCombinedDouble], type[KetCombinedComplexDouble]
    ]

    def __init__(
        self,
        systems: list[SystemAtomBase],
        m: Optional[tuple[float, float]] = None,
        energy: Union[tuple[float, float], tuple["PlainQuantity", "PlainQuantity"], None] = None,
        energy_unit: str = "pint",
    ) -> None:
        creator = self._BasisCombinedCreator()
        for system in systems:
            creator.add(system._cpp)
        if m is not None:
            creator.restrict_quantum_number_m(*m)
        if energy is not None:
            min_energy_au = Qty(energy[0], energy_unit).to_base("energy")
            max_energy_au = Qty(energy[1], energy_unit).to_base("energy")
            creator.restrict_energy(min_energy_au, max_energy_au)
        self._cpp = creator.create()

    @classmethod
    def _from_cpp_object(
        cls,
        cpp_obj: Union[
            _backend.BasisCombinedFloat,
            _backend.BasisCombinedComplexFloat,
            _backend.BasisCombinedDouble,
            _backend.BasisCombinedComplexDouble,
        ],
    ):
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    @cached_property
    def kets(
        self,
    ) -> list[Union[KetCombinedFloat, KetCombinedComplexFloat, KetCombinedDouble, KetCombinedComplexDouble]]:
        kets = [self._KetCombined._from_cpp_object(ket) for ket in self._cpp.get_kets()]
        return kets

    @property
    def number_of_states(self) -> int:
        return self._cpp.get_number_of_states()

    @property
    def number_of_kets(self) -> int:
        return self._cpp.get_number_of_kets()

    @property
    def coefficients(self) -> scipy.sparse.csr_matrix:
        return self._cpp.get_coefficients()


class BasisCombinedFloat(BasisCombinedBase):
    _cpp: _backend.BasisCombinedFloat
    _BasisCombinedCreator = _backend.BasisCombinedCreatorFloat
    _KetCombined = KetCombinedFloat


class BasisCombinedComplexFloat(BasisCombinedBase):
    _cpp: _backend.BasisCombinedComplexFloat
    _BasisCombinedCreator = _backend.BasisCombinedCreatorComplexFloat
    _KetCombined = KetCombinedComplexFloat


class BasisCombinedDouble(BasisCombinedBase):
    _cpp: _backend.BasisCombinedDouble
    _BasisCombinedCreator = _backend.BasisCombinedCreatorDouble
    _KetCombined = KetCombinedDouble


class BasisCombinedComplexDouble(BasisCombinedBase):
    _cpp: _backend.BasisCombinedComplexDouble
    _BasisCombinedCreator = _backend.BasisCombinedCreatorComplexDouble
    _KetCombined = KetCombinedComplexDouble
