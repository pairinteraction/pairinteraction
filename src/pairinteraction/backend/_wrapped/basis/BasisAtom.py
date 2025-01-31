from typing import TYPE_CHECKING, Any, ClassVar, Optional, TypeVar, Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.basis.Basis import BasisBase
from pairinteraction.backend._wrapped.database.Database import Database
from pairinteraction.backend._wrapped.get_functions import get_cpp_parity
from pairinteraction.backend._wrapped.ket.KetAtom import KetAtomBase, KetAtomDouble, KetAtomFloat
from pairinteraction.backend._wrapped.Parity import Parity
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

KetAtomType = TypeVar("KetAtomType", bound=KetAtomBase)
UnionCPPBasisAtom = Union[
    _backend.BasisAtomFloat, _backend.BasisAtomComplexFloat, _backend.BasisAtomDouble, _backend.BasisAtomComplexDouble
]
UnionTypeCPPBasisAtomCreator = Union[
    type[_backend.BasisAtomCreatorFloat],
    type[_backend.BasisAtomCreatorComplexFloat],
    type[_backend.BasisAtomCreatorDouble],
    type[_backend.BasisAtomCreatorComplexDouble],
]


class BasisAtomBase(BasisBase[KetAtomType]):
    _cpp: UnionCPPBasisAtom
    _cpp_creator: ClassVar[UnionTypeCPPBasisAtomCreator]

    def __init__(  # noqa: C901
        self,
        species: str,
        n: Optional[tuple[int, int]] = None,
        nu: Optional[tuple[float, float]] = None,
        l: Optional[tuple[float, float]] = None,
        s: Optional[tuple[float, float]] = None,
        j: Optional[tuple[float, float]] = None,
        f: Optional[tuple[float, float]] = None,
        m: Optional[tuple[float, float]] = None,
        energy: Union[tuple[float, float], tuple["PlainQuantity[float]", "PlainQuantity[float]"], None] = None,
        energy_unit: str = "pint",
        parity: Optional[Parity] = None,
        database: Optional[Database] = None,
        additional_kets: Optional[list[KetAtomType]] = None,
    ) -> None:
        """Create a basis for a single atom.

        Add all KetAtom objects that match the given quantum numbers to the basis.
        The initial coefficients matrix is a unit matrix, i.e. the first basis state is the first ket, etc.
        The BasisAtom coefficients matrix will always be square,
        i.e. the number of kets is equal to the number of states.

        Examples:
            >>> import pairinteraction.backend.double as pi
            >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
            >>> energy_min, energy_max = ket.get_energy(unit="GHz") - 100, ket.get_energy(unit="GHz") + 100
            >>> basis = pi.BasisAtom("Rb", n=(57, 63), l=(0, 3), energy=(energy_min, energy_max), energy_unit="GHz")
            >>> print(basis)
            BasisAtomDouble object with 140 states and 140 kets

        Args:
            species: The species of the atom.
            n: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            nu: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            l: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            s: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            j: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            f: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            m: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            energy: tuple of (min, max) value for the energy. Default None, i.e. add all available states.
            energy_unit: In which unit the energy values are given, e.g. "GHz".
                Default "pint", i.e. energy is provided as pint object.
            parity: The parity of the states to consider. Default None, i.e. add all available states.
            database: Which database to use. Default None, i.e. use the global database instance.
            additional_kets: List of additional kets to add to the basis. Default None.

        """
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
            min_energy_au = QuantityScalar(energy[0], energy_unit).to_base("ENERGY")
            max_energy_au = QuantityScalar(energy[1], energy_unit).to_base("ENERGY")
            creator.restrict_energy(min_energy_au, max_energy_au)
        if database is None:
            if Database.get_global_database() is None:
                Database.initialize_global_database()
            database = Database.get_global_database()
        if additional_kets is not None:
            for ket in additional_kets:
                creator.append_ket(ket._cpp)  # type: ignore [reportPrivateUsage]
        self._cpp = creator.create(database._cpp)  # type: ignore [reportPrivateUsage]


class BasisAtomFloat(BasisAtomBase[KetAtomFloat]):
    _cpp: _backend.BasisAtomFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisAtomCreatorFloat
    _TypeKet = KetAtomFloat


class BasisAtomComplexFloat(BasisAtomBase[KetAtomFloat]):
    _cpp: _backend.BasisAtomComplexFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisAtomCreatorComplexFloat
    _TypeKet = KetAtomFloat


class BasisAtomDouble(BasisAtomBase[KetAtomDouble]):
    _cpp: _backend.BasisAtomDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisAtomCreatorDouble
    _TypeKet = KetAtomDouble


class BasisAtomComplexDouble(BasisAtomBase[KetAtomDouble]):
    _cpp: _backend.BasisAtomComplexDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisAtomCreatorComplexDouble
    _TypeKet = KetAtomDouble


BasisAtom = BasisAtomBase[Any]
