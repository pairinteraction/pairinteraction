# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Literal, Optional, Union, overload

import numpy as np
from scipy.special import exprel

from pairinteraction import _backend
from pairinteraction._wrapped.database.database import Database
from pairinteraction._wrapped.enums import OperatorType, Parity, get_cpp_parity
from pairinteraction._wrapped.ket.ket import KetBase
from pairinteraction.units import QuantityArray, QuantityScalar, ureg

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.units import NDArray, PintArray, PintComplex, PintFloat


class KetAtom(KetBase):
    """Ket for an atomic basis state.

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
        >>> import pairinteraction.real as pi
        >>> ket_sqdt = pi.KetAtom("Rb", n=60, l=0, m=0.5)
        >>> (ket_sqdt.species, ket_sqdt.n, ket_sqdt.l, ket_sqdt.j, ket_sqdt.m, ket_sqdt.s)
        ('Rb', 60, 0.0, 0.5, 0.5, 0.5)
        >>> print(ket_sqdt)
        |Rb:60,S_1/2,1/2⟩
        >>> ket_mqdt = pi.KetAtom("Yb174_mqdt", nu=60, l=1, f=1, m=1)
        >>> (ket_mqdt.species, round(ket_mqdt.nu, 3), ket_mqdt.f, ket_mqdt.m)
        ('Yb174_mqdt', 60.049, 1.0, 1.0)
        >>> print(ket_mqdt)
        |Yb174:S=0.4,nu=60.0,L=1.0,J=1,1⟩

    """

    _cpp: _backend.KetAtom
    _cpp_creator = _backend.KetAtomCreator

    def __init__(  # noqa: C901, PLR0912
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
        energy: Union[float, "PintFloat", None] = None,
        energy_unit: Optional[str] = None,
        parity: Optional[Parity] = None,
        database: Optional[Database] = None,
    ) -> None:
        """Create a single atomic canonical basis state, which is defined by its species and quantum numbers.

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
            energy_au = QuantityScalar.from_pint_or_unit(energy, energy_unit, "energy").to_base_unit()
            creator.set_energy(energy_au)
        if f is not None:
            creator.set_quantum_number_f(f)
        if m is not None:
            creator.set_quantum_number_m(m)
        if parity is not None:
            creator.set_parity(get_cpp_parity(parity))
        if n is not None:
            if not (isinstance(n, int) or n.is_integer()):
                raise ValueError("Quantum number n must be an integer.")
            creator.set_quantum_number_n(int(n))
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
        self._cpp = creator.create(database._cpp)
        self._database = database

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, KetAtom):
            raise NotImplementedError
        if (self.species != other.species) or (self.m != other.m):
            return False
        return all(
            getattr(self, attr) == getattr(other, attr) for attr in ["n", "nu", "l", "s", "j", "l_ryd", "j_ryd", "f"]
        )

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
        """The expectation value of the total angular quantum number j_{Ryd} of the Rydberg electron."""
        return self._cpp.get_quantum_number_j_ryd()

    @property
    def nui_std(self) -> float:
        """The standard deviation of the effective principal quantum numbers nu_i of the channels."""
        return self._cpp.get_quantum_number_nui_std()

    @property
    def l_std(self) -> float:
        """The standard deviation of the orbital quantum number l of all valence electrons."""
        return self._cpp.get_quantum_number_l_std()

    @property
    def s_std(self) -> float:
        """The standard deviation of the total spin quantum number s of all valence electrons."""
        return self._cpp.get_quantum_number_s_std()

    @property
    def j_std(self) -> float:
        """The standard deviation of the total angular quantum number j of all valence electrons."""
        return self._cpp.get_quantum_number_j_std()

    @property
    def l_ryd_std(self) -> float:
        """The standard deviation of the orbital quantum number l_{Ryd} of the Rydberg electron."""
        return self._cpp.get_quantum_number_l_ryd_std()

    @property
    def j_ryd_std(self) -> float:
        """The standard deviation of the total angular quantum number j_{Ryd} of the Rydberg electron."""
        return self._cpp.get_quantum_number_j_ryd_std()

    @property
    def is_j_total_momentum(self) -> bool:
        """Whether j is the total momentum quantum number, otherwise f is the total momentum quantum number."""
        return self._cpp.is_j_total_momentum()

    @property
    def is_calculated_with_mqdt(self) -> bool:
        """Whether the state was calculated with multi-channel quantum defect theory."""
        return self._cpp.is_calculated_with_mqdt()

    @property
    def underspecified_channel_contribution(self) -> float:
        """The contribution of channels whose quantum numbers are not exactly known."""
        return self._cpp.get_underspecified_channel_contribution()

    @overload
    def get_matrix_element(
        self, ket: "Self", operator: OperatorType, q: int, unit: None = None
    ) -> Union["PintFloat", "PintComplex"]: ...  # type: ignore [type-var] # see "PintComplex"

    @overload
    def get_matrix_element(self, ket: "Self", operator: OperatorType, q: int, unit: str) -> Union[float, complex]: ...

    def get_matrix_element(
        self, ket: "Self", operator: OperatorType, q: int, unit: Optional[str] = None
    ) -> Union["PintFloat", "PintComplex", float, complex]:
        """Get the matrix element between two atomic basis states from the database.

        Args:
            ket: The second atomic basis state to calculate the matrix element with.
            operator: The operator, for which to calculate the matrix element.
            q: The index for the matrix element.
            unit: The unit to return the matrix element in. Default None will return a `pint.Quantity`.

        Returns:
            The matrix element between the two states in the given unit or as a `pint.Quantity`.

        """
        from pairinteraction._wrapped.basis.basis_atom import BasisAtomReal

        basis = BasisAtomReal(self.species, additional_kets=[self, ket], database=self.database)
        state_1 = basis.get_corresponding_state(self)
        state_2 = basis.get_corresponding_state(ket)

        return state_1.get_matrix_element(state_2, operator, q, unit=unit)

    @overload
    def get_spontaneous_transition_rates(self, unit: None = None) -> tuple[list["KetAtom"], "PintArray"]: ...

    @overload
    def get_spontaneous_transition_rates(self, unit: str) -> tuple[list["KetAtom"], "NDArray"]: ...

    def get_spontaneous_transition_rates(
        self, unit: Optional[str] = None
    ) -> tuple[list["KetAtom"], Union["NDArray", "PintArray"]]:
        """Calculate the spontaneous transition rates for the KetAtom.

        The spontaneous transition rates are given by the Einstein A coefficients.

        Args:
            unit: The unit to which to convert the result.
                Default None will return a `pint.Quantity`.

        Returns:
            The relevant states and the transition rates.

        """
        relevant_kets, transition_rates_au = self._get_transition_rates("spontaneous")
        transition_rates = QuantityArray.from_base_unit(transition_rates_au, "transition_rate").to_pint_or_unit(unit)
        return relevant_kets, transition_rates

    @overload
    def get_black_body_transition_rates(
        self, temperature: Union[float, "PintFloat"], temperature_unit: Optional[str] = None, unit: None = None
    ) -> tuple[list["KetAtom"], "PintArray"]: ...

    @overload
    def get_black_body_transition_rates(
        self, temperature: "PintFloat", *, unit: str
    ) -> tuple[list["KetAtom"], "NDArray"]: ...

    @overload
    def get_black_body_transition_rates(
        self, temperature: float, temperature_unit: str, unit: str
    ) -> tuple[list["KetAtom"], "NDArray"]: ...

    def get_black_body_transition_rates(
        self, temperature: Union[float, "PintFloat"], temperature_unit: Optional[str] = None, unit: Optional[str] = None
    ) -> tuple[list["KetAtom"], Union["NDArray", "PintArray"]]:
        """Calculate the black body transition rates of the KetAtom.

        The black body transitions rates are given by the Einstein B coefficients,
        with a weight factor given by Planck's law.

        Args:
            temperature: The temperature, for which to calculate the black body transition rates.
            temperature_unit: The unit of the temperature.
                Default None will assume the temperature is given as `pint.Quantity`.
            unit: The unit to which to convert the result.
                Default None will return a `pint.Quantity`.

        Returns:
            The relevant states and the transition rates.

        """
        temperature_au = QuantityScalar.from_pint_or_unit(temperature, temperature_unit, "temperature").to_base_unit()
        relevant_kets, transition_rates_au = self._get_transition_rates("black_body", temperature_au)
        transition_rates = QuantityArray.from_base_unit(transition_rates_au, "transition_rate").to_pint_or_unit(unit)
        return relevant_kets, transition_rates

    @overload
    def get_lifetime(
        self,
        temperature: Union[float, "PintFloat", None] = None,
        temperature_unit: Optional[str] = None,
        unit: None = None,
    ) -> "PintFloat": ...

    @overload
    def get_lifetime(self, *, unit: str) -> float: ...

    @overload
    def get_lifetime(self, temperature: "PintFloat", *, unit: str) -> float: ...

    @overload
    def get_lifetime(self, temperature: float, temperature_unit: str, unit: str) -> float: ...

    def get_lifetime(
        self,
        temperature: Union[float, "PintFloat", None] = None,
        temperature_unit: Optional[str] = None,
        unit: Optional[str] = None,
    ) -> Union[float, "PintFloat"]:
        """Calculate the lifetime of the KetAtom.

        The lifetime is the inverse of the sum of all transition rates.

        Args:
            temperature: The temperature, for which to calculate the black body transition rates.
                Default None will not include black body transitions.
            temperature_unit: The unit of the temperature.
                Default None will assume the temperature is given as `pint.Quantity`.
            unit: The unit to which to convert the result.
                Default None will return a `pint.Quantity`.

        Returns:
            The lifetime of the state.

        """
        _, transition_rates = self.get_spontaneous_transition_rates()
        transition_rates_au = transition_rates.to_base_units().magnitude
        if temperature is not None:
            _, black_body_transition_rates = self.get_black_body_transition_rates(temperature, temperature_unit)
            transition_rates_au = np.append(transition_rates_au, black_body_transition_rates.to_base_units().magnitude)

        lifetime_au = 1 / np.sum(transition_rates_au)

        return QuantityScalar.from_base_unit(lifetime_au, "time").to_pint_or_unit(unit)

    def _get_transition_rates(
        self, which_transitions: Literal["spontaneous", "black_body"], temperature_au: Union[float, None] = None
    ) -> tuple[list["KetAtom"], "NDArray"]:
        from pairinteraction._wrapped.basis.basis_atom import BasisAtomReal
        from pairinteraction._wrapped.system.system_atom import SystemAtomReal

        assert which_transitions in ["spontaneous", "black_body"]

        is_spontaneous = which_transitions == "spontaneous"
        n_max = self.n + 30

        energy_range = None
        if is_spontaneous:
            energy_range = (-1, self.get_energy("hartree"))

        basis = BasisAtomReal(
            self.species,
            n=(1, n_max),
            l=(self.l - 1, self.l + 1),
            m=(self.m - 1, self.m + 1),
            energy=energy_range,
            energy_unit="hartree",
            additional_kets=[self],  # needed to make get_matrix_elements(self, ...) work
            database=self.database,
        )
        system = SystemAtomReal(basis)

        relevant_kets = basis.kets
        energy_differences_au = self.get_energy("hartree") - system.get_eigenenergies("hartree")
        electric_dipole_moments_au = np.zeros(len(basis.kets))
        for q in [-1, 0, 1]:
            # the different entries are only at most once nonzero -> we can just add the arrays
            el_di_m = basis.get_matrix_elements(self, "electric_dipole", q)
            electric_dipole_moments_au += el_di_m.to_base_units().magnitude

        transition_rates_au = (
            (4 / 3)
            * np.abs(electric_dipole_moments_au) ** 2
            * energy_differences_au**2
            / ureg.Quantity(1, "speed_of_light").to_base_units().magnitude ** 3
        )

        if is_spontaneous:
            transition_rates_au *= energy_differences_au
        else:
            assert temperature_au is not None, "Temperature must be given for black body transitions."
            if temperature_au == 0:
                transition_rates_au *= 0
            else:  # for numerical stability we use 1 / exprel(x) = x / (exp(x) - 1)
                transition_rates_au *= temperature_au / exprel(energy_differences_au / temperature_au)

        mask = transition_rates_au != 0
        relevant_kets = [ket for ket, is_relevant in zip(relevant_kets, mask) if is_relevant]
        transition_rates_au = transition_rates_au[mask]
        return relevant_kets, transition_rates_au
