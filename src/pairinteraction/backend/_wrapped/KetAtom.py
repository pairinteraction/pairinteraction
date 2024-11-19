from typing import TYPE_CHECKING, Optional, Union

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.Ket import KetBase
from pairinteraction.unit_system import Qty

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

UnionCPPKetAtom = Union[_backend.KetAtomFloat, _backend.KetAtomDouble]
UnionTypeCPPKetAtomCreator = Union[type[_backend.KetAtomCreatorFloat], type[_backend.KetAtomCreatorDouble]]


class KetAtomBase(KetBase):
    _cpp: UnionCPPKetAtom
    _cpp_creator: UnionTypeCPPKetAtomCreator

    def __init__(
        self,
        species: str,
        n: Optional[int] = None,
        nu: Optional[float] = None,
        l: Optional[int] = None,
        s: Optional[float] = None,
        j: Optional[float] = None,
        f: Optional[float] = None,
        m: Optional[float] = None,
        energy: Union[float, "PlainQuantity", None] = None,
        energy_unit: str = "pint",
        parity: Optional[_backend.Parity] = None,
        database: Optional[_backend.Database] = None,
    ) -> None:
        creator = self._cpp_creator()
        creator.set_species(species)
        if n is not None:
            creator.set_quantum_number_n(n)
        if nu is not None:
            creator.set_quantum_number_nu(nu)
        if l is not None:
            creator.set_quantum_number_l(l)
        if s is not None:
            creator.set_quantum_number_s(s)
        if j is not None:
            creator.set_quantum_number_j(j)
        if f is not None:
            creator.set_quantum_number_f(f)
        if m is not None:
            creator.set_quantum_number_m(m)
        if energy is not None:
            energy_au = Qty(energy, energy_unit).to_base("energy")
            creator.set_energy(energy_au)
        if parity is not None:
            creator.set_parity(parity)
        if database is None:
            database = _backend.Database.get_global_instance(
                download_missing=False, wigner_in_memory=True, database_dir=""
            )  # Use the default database
        self._cpp = creator.create(database)

    @property
    def species(self) -> str:
        return self._cpp.get_species()

    @property
    def n(self) -> int:
        return self._cpp.get_quantum_number_n()

    @property
    def nu(self) -> float:
        return self._cpp.get_quantum_number_nu()

    @property
    def l(self) -> float:  # noqa: E743
        return self._cpp.get_quantum_number_l()

    @property
    def s(self) -> float:
        return self._cpp.get_quantum_number_s()

    @property
    def j(self) -> float:
        return self._cpp.get_quantum_number_j()


class KetAtomFloat(KetAtomBase):
    _cpp: _backend.KetAtomFloat
    _cpp_creator = _backend.KetAtomCreatorFloat


class KetAtomDouble(KetAtomBase):
    _cpp: _backend.KetAtomDouble
    _cpp_creator = _backend.KetAtomCreatorDouble
