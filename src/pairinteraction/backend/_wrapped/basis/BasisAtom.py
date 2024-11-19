from typing import TYPE_CHECKING, Optional, Union

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.basis.Basis import BasisBase
from pairinteraction.backend._wrapped.Database import Database
from pairinteraction.backend._wrapped.ket.KetAtom import KetAtomBase, KetAtomDouble, KetAtomFloat
from pairinteraction.backend._wrapped.Parity import Parity, get_cpp_parity
from pairinteraction.unit_system import Qty

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

UnionCPPBasisAtom = Union[
    _backend.BasisAtomFloat, _backend.BasisAtomComplexFloat, _backend.BasisAtomDouble, _backend.BasisAtomComplexDouble
]
UnionTypeCPPBasisAtomCreator = Union[
    type[_backend.BasisAtomCreatorFloat],
    type[_backend.BasisAtomCreatorComplexFloat],
    type[_backend.BasisAtomCreatorDouble],
    type[_backend.BasisAtomCreatorComplexDouble],
]
UnionTypeKetAtom = Union[type[KetAtomFloat], type[KetAtomDouble]]


class BasisAtomBase(BasisBase[KetAtomBase]):
    _cpp: UnionCPPBasisAtom
    _cpp_creator: UnionTypeCPPBasisAtomCreator
    _TypeKet: UnionTypeKetAtom

    def __init__(
        self,
        species: str,
        n: Optional[tuple[int, int]] = None,
        nu: Optional[tuple[float, float]] = None,
        l: Optional[tuple[float, float]] = None,
        s: Optional[tuple[float, float]] = None,
        j: Optional[tuple[float, float]] = None,
        f: Optional[tuple[float, float]] = None,
        m: Optional[tuple[float, float]] = None,
        energy: Union[tuple[float, float], tuple["PlainQuantity", "PlainQuantity"], None] = None,
        energy_unit: str = "pint",
        parity: Optional[Parity] = None,
        database: Optional[Database] = None,
    ) -> None:
        creator = self._cpp_creator()
        creator.set_species(species)
        if f is not None:
            creator.restrict_quantum_number_f(*f)
        if m is not None:
            creator.restrict_quantum_number_m(*m)
        if parity is not None:
            creator.restrict_parity(get_cpp_parity(parity))
        if n is not None:
            creator.restrict_quantum_number_n(*n)
        if nu is not None:
            creator.restrict_quantum_number_nu(*nu)
        if l is not None:
            creator.restrict_quantum_number_l(*l)
        if s is not None:
            creator.restrict_quantum_number_s(*s)
        if j is not None:
            creator.restrict_quantum_number_j(*j)
        if energy is not None:
            min_energy_au = Qty(energy[0], energy_unit).to_base("energy")
            max_energy_au = Qty(energy[1], energy_unit).to_base("energy")
            creator.restrict_energy(min_energy_au, max_energy_au)
        if database is None:
            database = Database.get_global_instance()
        self._cpp = creator.create(database._cpp)


class BasisAtomFloat(BasisAtomBase):
    _cpp: _backend.BasisAtomFloat
    _cpp_creator = _backend.BasisAtomCreatorFloat
    _TypeKet = KetAtomFloat


class BasisAtomComplexFloat(BasisAtomBase):
    _cpp: _backend.BasisAtomComplexFloat
    _cpp_creator = _backend.BasisAtomCreatorComplexFloat
    _TypeKet = KetAtomFloat


class BasisAtomDouble(BasisAtomBase):
    _cpp: _backend.BasisAtomDouble
    _cpp_creator = _backend.BasisAtomCreatorDouble
    _TypeKet = KetAtomDouble


class BasisAtomComplexDouble(BasisAtomBase):
    _cpp: _backend.BasisAtomComplexDouble
    _cpp_creator = _backend.BasisAtomCreatorComplexDouble
    _TypeKet = KetAtomDouble
