from typing import TYPE_CHECKING, ClassVar, Optional, Union, overload

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.database.Database import Database
from pairinteraction.backend._wrapped.get_functions import (
    get_basis_atom_class_from_ket,
    get_cpp_parity,
)
from pairinteraction.backend._wrapped.ket.Ket import KetBase
from pairinteraction.backend._wrapped.OperatorType import OperatorType
from pairinteraction.backend._wrapped.Parity import Parity
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity
    from typing_extensions import Self


UnionCPPKetAtom = Union[_backend.KetAtomFloat, _backend.KetAtomDouble]
UnionTypeCPPKetAtomCreator = Union[type[_backend.KetAtomCreatorFloat], type[_backend.KetAtomCreatorDouble]]


class KetAtomBase(KetBase):
    _cpp: UnionCPPKetAtom  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator: ClassVar[UnionTypeCPPKetAtomCreator]

    def __init__(  # noqa: C901
        self,
        species: str,
        n: Optional[int] = None,
        nu: Optional[float] = None,
        nui: Optional[float] = None,
        l: Optional[float] = None,
        s: Optional[float] = None,
        j: Optional[float] = None,
        l_ryd: Optional[float] = None,
        j_ryd: Optional[float] = None,
        f: Optional[float] = None,
        m: Optional[float] = None,
        energy: Union[float, "PlainQuantity[float]", None] = None,
        energy_unit: Optional[str] = None,
        parity: Optional[Parity] = None,
        database: Optional[Database] = None,
    ) -> None:
        """Create a single atomic canonical basis state, which is defined by its species and quantum numbers.

        Each KetAtom object uniquely represents a single atomic basis state
        (and therefore all KetAtom objects are orthogonal).
        When initializing a KetAtom you have to provide the species of the atom and a combination of quantum numbers,
        which uniquely define a single atomic basis state (this always includes providing a magnetic quantum number m).

        SQDT (Single Channel Quantum Defect Theory) for one valence electron (alkali atoms):
            The quantum numbers n (int), l (int), j (half-int) and m (half-int)
            should be used to define the desired atomic basis state.
            All other quantum numbers are trivially derived from these:
            s = 1/2, f = j (we neglect hyperfine interaction for SQDT),
            nu = n - delta, l_ryd = l, j_ryd = j.

        SQDT (Single Channel Quantum Defect Theory) for two valence electrons (earth-alkaline atoms):
            The quantum numbers n (int), l_ryd (int), j (int) and m (int)
            should be used to define the desired atomic basis state.
            The spin quantum number s is taken from the species label,
            which must end either with "_singlet" (s=0) or "_triplet" (s=1).
            Again we neglect hyperfine interaction, thus f = j. And nu = n - delta.
            All other quantum numbers are not necessarily eigenvalues anymore and are given as expectation values.

        MQDT (Multi Channel Quantum Defect Theory) for two valence electrons (earth-alkaline atoms):
            The quantum numbers nu (float), f (int or half-int) and m (int or half-int) are still good quantum numbers.
            All other quantum numbers (like l, s, j, l_ryd, j_ryd) are not necessarily eigenvalues anymore.
            You can still provide them to specify the atomic basis state,
            whose expectation value is closest to the provided value.

        Examples:
            >>> import pairinteraction.backend.double as pi
            >>> ket_sqdt = pi.KetAtom("Rb", n=60, l=0, m=0.5)
            >>> (ket_sqdt.species, ket_sqdt.n, ket_sqdt.l, ket_sqdt.j, ket_sqdt.m, ket_sqdt.s)
            ('Rb', 60, 0.0, 0.5, 0.5, 0.5)
            >>> print(ket_sqdt)
            Rb:60,S_1/2,1/2
            >>> ket_mqdt = pi.KetAtom("Yb174_mqdt", nu=60, l=1, f=1, m=1)
            >>> (ket_mqdt.species, round(ket_mqdt.nu, 3), ket_mqdt.f, ket_mqdt.m)
            ('Yb174_mqdt', 60.049, 1.0, 1.0)
            >>> print(ket_mqdt)
            Yb174:S=0.4,nu=60.0,L=1.0,J=1,1

        Args:
            species: See attribute.
            n: See attribute. Default None, i.e. load from the database.
            nu: See attribute. Default None, i.e. load from the database.
            nui: See attribute. Default None, i.e. load from the database.
            l: See attribute. Default None, i.e. load from the database.
            s: See attribute. Default None, i.e. load from the database.
            j: See attribute. Default None, i.e. load from the database.
            l_ryd: See attribute. Default None, i.e. load from the database.
            j_ryd: See attribute. Default None, i.e. load from the database.
            f: See attribute. Default None, i.e. load from the database.
            m: See attribute. This should always be provided.
            energy: See attribute. Default None, i.e. load from the database.
            energy_unit: In which unit the energy is given, e.g. "GHz".
                Default None, i.e. energy is provided as pint object.
            parity: See attribute. Default None, i.e. load from the database.
            database: Which database to use. Default None, i.e. use the global database instance.

        """
        creator = self._cpp_creator()
        creator.set_species(species)
        if energy is not None:
            energy_au = QuantityScalar(energy, energy_unit).to_base("ENERGY")
            creator.set_energy(energy_au)
        if f is not None:
            creator.set_quantum_number_f(f)
        if m is not None:
            creator.set_quantum_number_m(m)
        if parity is not None:
            creator.set_parity(get_cpp_parity(parity))
        if n is not None:
            creator.set_quantum_number_n(n)
        if nu is not None:
            creator.set_quantum_number_nu(nu)
        if nui is not None:
            creator.set_quantum_number_nui(nui)
        if l is not None:
            creator.set_quantum_number_l(l)
        if s is not None:
            creator.set_quantum_number_s(s)
        if j is not None:
            creator.set_quantum_number_j(j)
        if l_ryd is not None:
            creator.set_quantum_number_l_ryd(l_ryd)
        if j_ryd is not None:
            creator.set_quantum_number_j_ryd(j_ryd)
        if database is None:
            if Database.get_global_database() is None:
                Database.initialize_global_database()
            database = Database.get_global_database()
        self._cpp = creator.create(database._cpp)  # type: ignore [reportIncompatibleVariableOverride, reportPrivateUsage]
        self._database = database

    @property
    def database(self) -> Database:
        """The database used for this object."""
        return self._database

    @property
    def species(self) -> str:
        """The atomic species."""
        return self._cpp.get_species()

    @property
    def n(self) -> int:
        """The principal quantum number n."""
        return self._cpp.get_quantum_number_n()

    @property
    def nu(self) -> float:
        """The effective principal quantum number nu."""
        return self._cpp.get_quantum_number_nu()

    @property
    def nui(self) -> float:
        """The expectation value of the effective principal quantum numbers nu_i of the channels."""
        return self._cpp.get_quantum_number_nui()

    @property
    def l(self) -> float:  # noqa: E743
        """The expectation value of the orbital quantum number l of all valence electrons."""
        return self._cpp.get_quantum_number_l()

    @property
    def s(self) -> float:
        """The expectation value of the total spin quantum number s of all valence electrons."""
        return self._cpp.get_quantum_number_s()

    @property
    def j(self) -> float:
        """The expectation value of the total angular quantum number j of all valence electrons."""
        return self._cpp.get_quantum_number_j()

    @property
    def l_ryd(self) -> float:
        """The expectation value of the orbital quantum number l_{Ryd} of the Rydberg electron."""
        return self._cpp.get_quantum_number_l_ryd()

    @property
    def j_ryd(self) -> float:
        """The expectation value of the total angular quantum number j_{Ryd} the Rydberg electron."""
        return self._cpp.get_quantum_number_j_ryd()

    @overload
    def get_matrix_element(self, ket: "Self", operator: OperatorType, q: int) -> "PlainQuantity[float]": ...

    @overload
    def get_matrix_element(self, ket: "Self", operator: OperatorType, q: int, unit: str) -> "float": ...

    def get_matrix_element(
        self,
        ket: "Self",
        operator: OperatorType,
        q: int,
        unit: Optional[str] = None,
    ):
        BasisAtomClass = get_basis_atom_class_from_ket(self)
        basis = BasisAtomClass(self.species, additional_kets=[self, ket], database=self.database)
        state_1 = basis.get_corresponding_state(self)

        matrixelements = state_1.get_matrix_elements(ket, operator, q, unit=unit)
        return matrixelements[0]


class KetAtomFloat(KetAtomBase):
    _cpp: _backend.KetAtomFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.KetAtomCreatorFloat


class KetAtomDouble(KetAtomBase):
    _cpp: _backend.KetAtomDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.KetAtomCreatorDouble


KetAtom = KetAtomBase
