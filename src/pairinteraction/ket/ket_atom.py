# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from functools import cached_property
from typing import TYPE_CHECKING, Literal, overload

import numpy as np
from scipy.special import exprel

from pairinteraction import _backend
from pairinteraction.database import Database
from pairinteraction.enums import OperatorType, Parity, int_to_parity, parity_to_int
from pairinteraction.ket.ket_base import KetBase
from pairinteraction.ket.utils import format_half_integer, get_l_label
from pairinteraction.units import QuantityArray, QuantityScalar, ureg

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.enums import OperatorType, Parity
    from pairinteraction.state import StateAtom, StateAtomReal
    from pairinteraction.units import NDArray, PintArray, PintComplex, PintFloat


logger = logging.getLogger(__name__)


class KetAtom(KetBase):
    """Ket for an atomic basis state.

    Each KetAtom object uniquely represents a single-atom basis state
    (and therefore all KetAtom objects are orthogonal).
    When initializing a KetAtom you have to provide the species of the atom and a combination of quantum numbers,
    which uniquely define a single-atom basis state (this always includes providing a magnetic quantum number m).

    SQDT (Single Channel Quantum Defect Theory) for one valence electron (alkali atoms):
        The quantum numbers n (int), l (int), j (half-int) and m (half-int)
        should be used to define the desired atomic basis state.
        All other quantum numbers are trivially derived from these:
        s = 1/2, f = j (we neglect hyperfine interaction for SQDT),
        nu = n - delta, l_ryd = l, j_ryd = j.

    SQDT (Single Channel Quantum Defect Theory) for two valence electrons (divalent atoms):
        The quantum numbers n (int), l_ryd (int), s (0 or 1), j (int) and m (int)
        should be used to define the desired atomic basis state.
        For divalent atoms the spin quantum number s selects the singlet (s=0) or the triplet (s=1) sector,
        which are both contained in the same database (e.g. "Sr88_sqdt").
        Again we neglect hyperfine interaction, thus f = j. And nu = n - delta.
        All other quantum numbers are not necessarily eigenvalues anymore and are given as expectation values.

    MQDT (Multi Channel Quantum Defect Theory) for two valence electrons (divalent atoms):
        The quantum numbers nu (float), f (int or half-int) and m (int or half-int) are still good quantum numbers.
        All other quantum numbers (like l, s, j, l_ryd, j_ryd) are not necessarily eigenvalues anymore.
        You can still provide them to specify the atomic basis state,
        whose expectation value is closest to the provided value.

    Examples:
        >>> import pairinteraction as pi
        >>> ket_s = pi.KetAtom("Rb", n=60, l=0, m=0.5)
        >>> (ket_s.species, ket_s.n, ket_s.l, ket_s.j, ket_s.m, ket_s.s)
        ('Rb', 60, 0.0, 0.5, 0.5, 0.5)
        >>> print(ket_s)
        |Rb:60,S_1/2,1/2⟩
        >>> print(ket_s.to_state())
        1.00 |Rb:60,S_1/2,1/2⟩
        >>> ket_p = pi.KetAtom("Rb", n=60, l=1, j=0.5, m=0.5)
        >>> print((2 * ket_p - ket_s).normalize())
        0.89 |Rb:60,P_1/2,1/2⟩ - 0.45 |Rb:60,S_1/2,1/2⟩
        >>> ket_mqdt = pi.KetAtom("Yb174_mqdt", nu=60, l=1, f=1, m=1)
        >>> (ket_mqdt.species, round(ket_mqdt.nu, 3), ket_mqdt.f, ket_mqdt.m)
        ('Yb174_mqdt', 60.049, 1.0, 1.0)
        >>> print(ket_mqdt)
        |Yb174:S=0.0,nu=60.0,L=1.0,J=1,1⟩

    """

    _cpp: _backend.KetAtom

    def __init__(
        self,
        species: str,
        n: int | None = None,
        nu: float | None = None,
        nui: float | None = None,
        l: float | None = None,
        s: float | None = None,
        j: float | None = None,
        l_ryd: float | None = None,
        j_ryd: float | None = None,
        f: float | None = None,
        m: float | None = None,
        energy: float | PintFloat | None = None,
        energy_unit: str | None = None,
        parity: Parity | None = None,
        database: Database | None = None,
    ) -> None:
        """Create a single-atom canonical basis state, which is defined by its species and quantum numbers.

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
        creator = _backend.KetAtomCreator()
        creator.set_species(species)
        if energy is not None:
            energy_au = QuantityScalar.convert_user_to_au(energy, energy_unit, "energy")
            creator.set_energy(energy_au)
        if n is not None and not (isinstance(n, int) or n.is_integer()):
            raise ValueError("Quantum number n must be an integer.")
        quantum_numbers = {
            "f": f,
            "m": m,
            "n": n,
            "nu": nu,
            "nui": nui,
            "l": l,
            "s": s,
            "j": j,
            "l_ryd": l_ryd,
            "j_ryd": j_ryd,
            "parity": parity_to_int(parity) if parity is not None else None,
        }
        for name, value in quantum_numbers.items():
            if value is not None:
                creator.set_quantum_number(name, value)
        if database is None:
            if Database.get_global_database() is None:
                Database.initialize_global_database()
            database = Database.get_global_database()
        try:
            self._cpp = creator.create(database._cpp)
        except _backend.KetNotUniqueError as err:
            candidates = [type(self)._from_cpp_object(ket) for ket in err.kets]  # type: ignore [attr-defined]
            labels = "\n".join(ket.get_label("ket") for ket in candidates)
            raise ValueError(f"The ket is not uniquely specified. Possible kets are:\n{labels}") from None

    def _get_raw_label(self) -> str:
        s, l, f, m = self.s, self.l, self.f, self.m

        label = self.species.split("_", 1)[0]
        label = label[0].upper() + label[1:]

        if not self.species.endswith("_mqdt"):
            if s == 0:
                label += "_singlet"
            elif s == 1:
                label += "_triplet"
            elif s != 0.5:
                logger.error("Unexpected spin quantum number s=%f for species %s.", s, self.species)

        label += ":"

        if self.species.endswith("_mqdt"):
            label += f"S={s:.1f},nu={self.nu:.1f},L={l:.1f},"
            label += "J=" if self.is_j_total_momentum else "F="
        else:
            label += f"{self.n:d},"
            label += get_l_label(l)
            label += "_"

        label += format_half_integer(f)
        label += "," + format_half_integer(m)

        return label

    @cached_property
    def database(self) -> Database:
        """The database from which the KetAtom was loaded."""
        database_cpp = self._cpp.get_database()
        return Database._from_cpp_object(database_cpp)

    @property
    def m(self) -> float:
        """The magnetic quantum number m (int or half-int)."""
        return self._cpp.get_quantum_number("m")

    @property
    def f(self) -> float:
        """The total momentum quantum number f (int or half-int)."""
        return self._cpp.get_quantum_number("f")

    @property
    def parity(self) -> Parity:
        """The parity of the ket."""
        return int_to_parity(int(self._cpp.get_quantum_number("parity")))

    @property
    def species(self) -> str:
        """The atomic species."""
        return self._cpp.get_species()

    @property
    def n(self) -> int:
        """The principal quantum number n."""
        return int(self._cpp.get_quantum_number("n"))

    @property
    def nu(self) -> float:
        """The effective principal quantum number nu."""
        return self._cpp.get_quantum_number("nu")

    @property
    def nui(self) -> float:
        """The expectation value of the effective principal quantum numbers nu_i of the channels."""
        return self._cpp.get_quantum_number("nui")

    @property
    def l(self) -> float:  # noqa: E743
        """The expectation value of the orbital quantum number l of all valence electrons."""
        return self._cpp.get_quantum_number("l")

    @property
    def s(self) -> float:
        """The expectation value of the total spin quantum number s of all valence electrons."""
        return self._cpp.get_quantum_number("s")

    @property
    def j(self) -> float:
        """The expectation value of the total angular quantum number j of all valence electrons."""
        return self._cpp.get_quantum_number("j")

    @property
    def l_ryd(self) -> float:
        """The expectation value of the orbital quantum number l_{Ryd} of the Rydberg electron."""
        return self._cpp.get_quantum_number("l_ryd")

    @property
    def j_ryd(self) -> float:
        """The expectation value of the total angular quantum number j_{Ryd} of the Rydberg electron."""
        return self._cpp.get_quantum_number("j_ryd")

    @property
    def nui_std(self) -> float:
        """The standard deviation of the effective principal quantum numbers nu_i of the channels."""
        return self._cpp.get_quantum_number_std("nui")

    @property
    def l_std(self) -> float:
        """The standard deviation of the orbital quantum number l of all valence electrons."""
        return self._cpp.get_quantum_number_std("l")

    @property
    def s_std(self) -> float:
        """The standard deviation of the total spin quantum number s of all valence electrons."""
        return self._cpp.get_quantum_number_std("s")

    @property
    def j_std(self) -> float:
        """The standard deviation of the total angular quantum number j of all valence electrons."""
        return self._cpp.get_quantum_number_std("j")

    @property
    def l_ryd_std(self) -> float:
        """The standard deviation of the orbital quantum number l_{Ryd} of the Rydberg electron."""
        return self._cpp.get_quantum_number_std("l_ryd")

    @property
    def j_ryd_std(self) -> float:
        """The standard deviation of the total angular quantum number j_{Ryd} of the Rydberg electron."""
        return self._cpp.get_quantum_number_std("j_ryd")

    @property
    def is_j_total_momentum(self) -> bool:
        """Whether j is the total momentum quantum number, otherwise f is the total momentum quantum number."""
        return bool(self._cpp.get_quantum_number("is_j_total_momentum"))

    @property
    def is_calculated_with_mqdt(self) -> bool:
        """Whether the state was calculated with multi-channel quantum defect theory."""
        return bool(self._cpp.get_quantum_number("is_calculated_with_mqdt"))

    @property
    def underspecified_channel_contribution(self) -> float:
        """The contribution of channels whose quantum numbers are not exactly known."""
        return self._cpp.get_quantum_number("underspecified_channel_contribution")

    def to_state(self) -> StateAtom:
        """Create a canonical state representing the single ket.

        The returned state has a minimal basis consisting only of this ket and a single coefficient equal to one.

        Returns:
            A state object representing the ket.

        """
        from pairinteraction.state import StateAtom

        return StateAtom([1], [self])

    def __add__(self, other: KetAtom | StateAtom) -> StateAtom:
        """Build the superposition of this ket and another ket or state.

        The ket is converted to a state via :meth:`to_state`, the resulting superposition is in general not normalized.
        """
        return self.to_state() + other

    def __sub__(self, other: KetAtom | StateAtom) -> StateAtom:
        """Build the superposition of this ket and the negative of another ket or state.

        The ket is converted to a state via :meth:`to_state`, the resulting superposition is in general not normalized.
        """
        return self.to_state() - other

    def __mul__(self, factor: complex) -> StateAtom:
        """Scale the ket by a complex amplitude, e.g. to build superpositions like ``2 * ket_s - 1j * ket_p``."""
        return self.to_state() * factor

    def __truediv__(self, factor: complex) -> StateAtom:
        """Scale the ket by the inverse of a complex amplitude."""
        return self.to_state() / factor

    def __neg__(self) -> StateAtom:
        """Flip the sign of the amplitude of the ket."""
        return -self.to_state()

    __rmul__ = __mul__  # for reverse multiplication, i.e. scalar * ket will use ket.__rmul__

    @overload
    def get_matrix_element(
        self, ket: Self | StateAtom, operator: OperatorType, q: int, unit: None = None
    ) -> PintFloat | PintComplex: ...

    @overload
    def get_matrix_element(
        self, ket: Self | StateAtom, operator: OperatorType, q: int, unit: str
    ) -> float | complex: ...

    def get_matrix_element(
        self, ket: Self | StateAtom, operator: OperatorType, q: int, unit: str | None = None
    ) -> PintFloat | PintComplex | float | complex:
        """Get the matrix element between two atomic basis states from the database.

        Args:
            ket: The second atomic basis state to calculate the matrix element with.
            operator: The operator, for which to calculate the matrix element.
            q: The index for the matrix element.
            unit: The unit to return the matrix element in. Default None will return a `pint.Quantity`.

        Returns:
            The matrix element between the two states in the given unit or as a `pint.Quantity`.

        """
        return self.to_state().get_matrix_element(ket, operator, q, unit=unit)

    @overload
    def get_spontaneous_transition_rates(self, unit: None = None) -> tuple[list[KetAtom], PintArray]: ...

    @overload
    def get_spontaneous_transition_rates(self, unit: str) -> tuple[list[KetAtom], NDArray]: ...

    def get_spontaneous_transition_rates(self, unit: str | None = None) -> tuple[list[KetAtom], NDArray | PintArray]:
        """Calculate the spontaneous transition rates for the KetAtom.

        The spontaneous transition rates are given by the Einstein A coefficients.

        Args:
            unit: The unit to which to convert the result.
                Default None will return a `pint.Quantity`.

        Returns:
            The relevant states and the transition rates.

        """
        relevant_kets, transition_rates_au = self._get_transition_rates("spontaneous")
        transition_rates = QuantityArray.convert_au_to_user(transition_rates_au, "transition_rate", unit)
        return relevant_kets, transition_rates

    @overload
    def get_black_body_transition_rates(
        self, temperature: float | PintFloat, temperature_unit: str | None = None, unit: None = None
    ) -> tuple[list[KetAtom], PintArray]: ...

    @overload
    def get_black_body_transition_rates(
        self, temperature: PintFloat, *, unit: str
    ) -> tuple[list[KetAtom], NDArray]: ...

    @overload
    def get_black_body_transition_rates(
        self, temperature: float, temperature_unit: str, unit: str
    ) -> tuple[list[KetAtom], NDArray]: ...

    def get_black_body_transition_rates(
        self, temperature: float | PintFloat, temperature_unit: str | None = None, unit: str | None = None
    ) -> tuple[list[KetAtom], NDArray | PintArray]:
        """Calculate the black body transition rates of the KetAtom.

        The black body transition rates are given by the Einstein B coefficients,
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
        temperature_au = QuantityScalar.convert_user_to_au(temperature, temperature_unit, "temperature")
        relevant_kets, transition_rates_au = self._get_transition_rates("black_body", temperature_au)
        transition_rates = QuantityArray.convert_au_to_user(transition_rates_au, "transition_rate", unit)
        return relevant_kets, transition_rates

    @overload
    def get_lifetime(
        self,
        temperature: float | PintFloat | None = None,
        temperature_unit: str | None = None,
        unit: None = None,
    ) -> PintFloat: ...

    @overload
    def get_lifetime(self, *, unit: str) -> float: ...

    @overload
    def get_lifetime(self, temperature: PintFloat, *, unit: str) -> float: ...

    @overload
    def get_lifetime(self, temperature: float, temperature_unit: str, unit: str) -> float: ...

    def get_lifetime(
        self,
        temperature: float | PintFloat | None = None,
        temperature_unit: str | None = None,
        unit: str | None = None,
    ) -> float | PintFloat:
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

        return QuantityScalar.convert_au_to_user(lifetime_au, "time", unit)

    def _get_transition_rates(
        self, which_transitions: Literal["spontaneous", "black_body"], temperature_au: float | None = None
    ) -> tuple[list[KetAtom], NDArray]:
        if not isinstance(self, KetAtomReal):
            from pairinteraction.basis import BasisAtom

            basis_atom_class = BasisAtom
        else:
            from pairinteraction.basis import BasisAtomReal

            basis_atom_class = BasisAtomReal

        assert which_transitions in ["spontaneous", "black_body"]
        is_spontaneous = which_transitions == "spontaneous"
        n_max = self.n + 30

        energy_range = None
        if is_spontaneous:
            energy_range = (-1, self.get_energy("hartree"))

        basis = basis_atom_class(
            self.species,
            n=(1, n_max),
            l=(self.l - 1, self.l + 1),
            m=(self.m - 1, self.m + 1),
            energy=energy_range,
            energy_unit="hartree",
            additional_kets=[self],  # needed to make get_matrix_elements(self, ...) work
            database=self.database,
        )

        relevant_kets = basis.kets
        energy_differences_au = np.abs(
            self.get_energy("hartree") - np.array([ket_cpp.get_energy() for ket_cpp in basis._cpp.get_kets()])
        )
        electric_dipole_moments_au = np.zeros(len(basis.kets), dtype=complex)
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
        relevant_kets = [ket for ket, is_relevant in zip(relevant_kets, mask, strict=True) if is_relevant]
        transition_rates_au = transition_rates_au[mask]
        return relevant_kets, transition_rates_au


class KetAtomReal(KetAtom):
    def to_state(self) -> StateAtomReal:
        from pairinteraction.state import StateAtomReal

        return StateAtomReal([1], [self])
