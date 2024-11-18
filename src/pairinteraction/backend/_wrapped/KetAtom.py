from typing import TYPE_CHECKING, Optional, Union

import pairinteraction.backend._backend as _backend
from pairinteraction.unit_system import Qty

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity


class KetAtomBase:
    _KetAtomCreator: Union[type[_backend.KetAtomCreatorFloat], type[_backend.KetAtomCreatorDouble]]

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
        creator = self._KetAtomCreator()
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

    @classmethod
    def _from_cpp_object(cls, cpp_obj: Union[_backend.KetAtomFloat, _backend.KetAtomDouble]):
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    def __str__(self) -> str:
        return self.label

    # # # Define convenience properties # # #
    @property
    def label(self) -> str:
        return self._cpp.get_label()

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

    @property
    def f(self) -> float:
        return self._cpp.get_quantum_number_f()

    @property
    def m(self) -> float:
        return self._cpp.get_quantum_number_m()

    @property
    def species(self) -> str:
        return self._cpp.get_species()

    def get_energy(self, unit: str = "pint"):
        energy_au = self._cpp.get_energy()
        energy = Qty.from_base(energy_au, "energy")
        return energy.to_unit(unit)


class KetAtomFloat(KetAtomBase):
    _cpp: _backend.KetAtomFloat
    _KetAtomCreator = _backend.KetAtomCreatorFloat


class KetAtomDouble(KetAtomBase):
    _cpp: _backend.KetAtomDouble
    _KetAtomCreator = _backend.KetAtomCreatorDouble
