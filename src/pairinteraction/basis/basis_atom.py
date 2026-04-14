# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from typing import TYPE_CHECKING, Any, overload

import numpy as np

from pairinteraction import _backend
from pairinteraction.basis.basis_base import BasisBase
from pairinteraction.database import Database
from pairinteraction.enums import get_cpp_operator_type, get_cpp_parity
from pairinteraction.ket import KetAtom
from pairinteraction.state import StateAtom, StateAtomReal
from pairinteraction.units import QuantityArray, QuantityScalar, QuantitySparse

if TYPE_CHECKING:
    from collections.abc import Sequence

    from scipy.sparse import csr_matrix
    from typing_extensions import Self

    from pairinteraction.enums import OperatorType, Parity
    from pairinteraction.units import NDArray, PintArray, PintFloat, PintSparse


class BasisAtom(BasisBase[KetAtom, StateAtom]):
    """Basis for a single atom.

    Add all KetAtom objects that match the given quantum numbers to the basis.
    The initial coefficients matrix is a unit matrix, i.e. the first basis state is the first ket, etc.
    The BasisAtom coefficients matrix will always be square,
    i.e. the number of kets is equal to the number of states.

    Examples:
        >>> import pairinteraction as pi
        >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
        >>> energy_min, energy_max = ket.get_energy(unit="GHz") - 100, ket.get_energy(unit="GHz") + 100
        >>> basis = pi.BasisAtom("Rb", n=(57, 63), l=(0, 3), energy=(energy_min, energy_max), energy_unit="GHz")
        >>> print(basis)
        BasisAtom(n=(57, 63), l=(0, 3), energy=(1008911.9216, 1009111.9216), energy_unit=GHz)

    """

    _cpp: _backend.BasisAtomComplex
    _cpp_creator = _backend.BasisAtomCreatorComplex
    _ket_class = KetAtom
    _state_class = StateAtom

    _parameters_creator: dict[str, Any] | None = None

    def __init__(  # noqa: C901, PLR0912, PLR0915
        self,
        species: str,
        n: tuple[int, int] | None = None,
        nu: tuple[float, float] | None = None,
        nui: tuple[float, float] | None = None,
        l: tuple[float, float] | None = None,
        s: tuple[float, float] | None = None,
        j: tuple[float, float] | None = None,
        l_ryd: tuple[float, float] | None = None,
        j_ryd: tuple[float, float] | None = None,
        f: tuple[float, float] | None = None,
        m: tuple[float, float] | None = None,
        energy: tuple[float, float] | tuple[PintFloat, PintFloat] | None = None,
        energy_unit: str | None = None,
        parity: Parity | None = None,
        database: Database | None = None,
        additional_kets: Sequence[KetAtom] | None = None,
    ) -> None:
        """Create a basis for a single atom.

        Args:
            species: The species of the atom.
            n: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            nu: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            nui: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            l: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            s: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            j: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            l_ryd: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            j_ryd: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            f: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            m: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            energy: tuple of (min, max) value for the energy. Default None, i.e. add all available states.
            energy_unit: In which unit the energy values are given, e.g. "GHz".
                Default None, i.e. energy is provided as pint object.
            parity: The parity of the states to consider. Default None, i.e. add all available states.
            database: Which database to use. Default None, i.e. use the global database instance.
            additional_kets: List of additional kets to add to the basis. Default None.

        """
        super().__init__()
        creator = self._cpp_creator()
        creator.set_species(species)

        self._parameters_creator = {"qns": {}}
        if n is not None:
            if not all(isinstance(x, int) or x.is_integer() for x in n):
                raise ValueError("Quantum numbers n must be integers.")
            n = (int(n[0]), int(n[1]))
            creator.restrict_quantum_number_n(*n)
            self._parameters_creator["qns"]["n"] = n
        if nu is not None:
            creator.restrict_quantum_number_nu(*nu)
            self._parameters_creator["qns"]["nu"] = nu
        if nui is not None:
            creator.restrict_quantum_number_nui(*nui)
            self._parameters_creator["qns"]["nui"] = nui
        if s is not None:
            creator.restrict_quantum_number_s(*s)
            self._parameters_creator["qns"]["s"] = s
        if l is not None:
            creator.restrict_quantum_number_l(*l)
            self._parameters_creator["qns"]["l"] = l
        if j is not None:
            self._parameters_creator["qns"]["j"] = j
            creator.restrict_quantum_number_j(*j)
        if l_ryd is not None:
            creator.restrict_quantum_number_l_ryd(*l_ryd)
            self._parameters_creator["qns"]["l_ryd"] = l_ryd
        if j_ryd is not None:
            self._parameters_creator["qns"]["j_ryd"] = j_ryd
            creator.restrict_quantum_number_j_ryd(*j_ryd)
        if f is not None:
            creator.restrict_quantum_number_f(*f)
            self._parameters_creator["qns"]["f"] = f
        if m is not None:
            creator.restrict_quantum_number_m(*m)
            self._parameters_creator["qns"]["m"] = m
        if parity is not None:
            creator.restrict_parity(get_cpp_parity(parity))
        if energy is not None:
            min_energy_au = QuantityScalar.convert_user_to_au(energy[0], energy_unit, "energy")
            max_energy_au = QuantityScalar.convert_user_to_au(energy[1], energy_unit, "energy")
            creator.restrict_energy(min_energy_au, max_energy_au)
            self._parameters_creator["energy"] = energy
            self._parameters_creator["energy_unit"] = energy_unit
        if database is None:
            if Database.get_global_database() is None:
                Database.initialize_global_database()
            database = Database.get_global_database()
        if additional_kets is not None:
            for ket in additional_kets:
                creator.append_ket(ket._cpp)
            self._parameters_creator["additional_kets"] = additional_kets
        self._cpp = creator.create(database._cpp)

    @classmethod
    def from_kets(
        cls: type[Self],
        kets: KetAtom | Sequence[KetAtom],
        delta_n: int | None = None,
        delta_nu: float | None = None,
        delta_nui: float | None = None,
        delta_l: float | None = None,
        delta_s: float | None = None,
        delta_j: float | None = None,
        delta_l_ryd: float | None = None,
        delta_j_ryd: float | None = None,
        delta_f: int | None = None,
        delta_m: int | None = None,
        delta_energy: float | PintFloat | None = None,
        delta_energy_unit: str | None = None,
        parity: Parity | None = None,
        database: Database | None = None,
        additional_kets: Sequence[KetAtom] | None = None,
    ) -> Self:
        """Create a BasisAtom from one or more kets and quantum number deltas.

        Currently a single big basis including all kets for the quantum numbers
        from min_value - delta to max_value + delta is returned.
        In the future this might change to return a basis including all states around the given kets +/- delta,
        but not necessarily all states between the given kets.

        For each quantum number, pass the corresponding ``delta_*`` argument to include
        all states within ``[min_value - delta, max_value + delta]``, where
        ``min_value`` / ``max_value`` are the extremes across all provided kets.
        If no ``delta_*`` is given for a quantum number, that quantum number is left
        unrestricted.

        Args:
            kets: The ket(s) around which the basis should be centered.
            delta_n: Half-width of the n window (integer steps).
                Default None means no restriction on n.
            delta_nu: Half-width of the nu window.
                Default None means no restriction on nu.
            delta_nui: Half-width of the nui window.
                Default None means no restriction on nui.
            delta_l: Half-width of the l window.
                Default None means no restriction on l.
            delta_s: Half-width of the s window.
                Default None means no restriction on s.
            delta_j: Half-width of the j window.
                Default None means no restriction on j.
            delta_l_ryd: Half-width of the l_ryd window.
                Default None means no restriction on l_ryd.
            delta_j_ryd: Half-width of the j_ryd window.
                Default None means no restriction on j_ryd.
            delta_f: Half-width of the f window (integer steps).
                Default None means no restriction on f.
            delta_m: Half-width of the m window (integer steps).
                Default None means no restriction on m.
            delta_energy: Half-width of the energy window around the energies of the
                provided kets. Default None means no energy restriction.
            delta_energy_unit: Unit for ``delta_energy`` (e.g. ``"GHz"``).
                Default None means pint quantities are used.
            parity: Restrict to states with this parity.
                Default None means no parity restriction.
            database: Database instance to use.
                Default None uses the global database.
            additional_kets: Extra kets to force-include in the basis.
                Default None.

        Returns:
            A new :class:`BasisAtom` centered around the provided kets.

        Examples:
            >>> import pairinteraction as pi
            >>> ket1 = pi.KetAtom("Rb", n=60, l=0, m=0.5)
            >>> ket2 = pi.KetAtom("Rb", n=59, l=0, m=0.5)
            >>> basis = pi.BasisAtom.from_kets([ket1, ket2], delta_n=2, delta_l=1)
            >>> basis.species
            'Rb'
            >>> all(57 <= k.n <= 62 for k in basis.kets)
            True

        """
        if isinstance(kets, KetAtom):
            kets = [kets]
        kets = list(kets)
        if len(kets) == 0:
            raise ValueError("kets must not be empty.")
        if len({ket.species for ket in kets}) > 1:
            raise ValueError(f"All kets must have the same species, but got: {sorted({ket.species for ket in kets})}.")

        def get_range(name: str, delta: float | None) -> tuple[float, float] | None:
            if delta is None:
                return None
            if name == "energy":
                values = [ket.get_energy(unit=delta_energy_unit) for ket in kets]
            else:
                values = [getattr(ket, name) for ket in kets]
            return (min(values) - delta, max(values) + delta)

        return cls(
            species=kets[0].species,
            n=get_range("n", delta_n),  # type: ignore [arg-type]
            nu=get_range("nu", delta_nu),
            nui=get_range("nui", delta_nui),
            l=get_range("l", delta_l),
            s=get_range("s", delta_s),
            j=get_range("j", delta_j),
            l_ryd=get_range("l_ryd", delta_l_ryd),
            j_ryd=get_range("j_ryd", delta_j_ryd),
            f=get_range("f", delta_f),
            m=get_range("m", delta_m),
            energy=get_range("energy", delta_energy),  # type: ignore [arg-type]
            energy_unit=delta_energy_unit,
            parity=parity,
            database=database,
            additional_kets=additional_kets,
        )

    def __repr__(self) -> str:
        if self._parameters_creator is None:
            return super().__repr__()

        args = ""
        if qns := self._parameters_creator.get("qns", None):
            args += ", ".join(f"{k}={v}" for k, v in qns.items())
        if (energy := self._parameters_creator.get("energy", None)) and (
            energy_unit := self._parameters_creator.get("energy_unit", None)
        ):
            args += f", energy=({energy[0]:.4f}, {energy[1]:.4f}), energy_unit={energy_unit}"
        if additional_kets := self._parameters_creator.get("additional_kets", None):
            args += f", additional_kets={additional_kets}"
        return f"{type(self).__name__}({args})"

    @property
    def database(self) -> Database:
        """The database used for this object."""
        return self.get_ket(0).database

    @property
    def species(self) -> str:
        """The atomic species."""
        return self.get_ket(0).species

    @overload
    def get_amplitudes(self, other: KetAtom | StateAtom) -> NDArray: ...

    @overload
    def get_amplitudes(self, other: Self) -> csr_matrix: ...

    def get_amplitudes(self, other: KetAtom | StateAtom | Self) -> NDArray | csr_matrix:
        if isinstance(other, KetAtom):
            return np.array(self._cpp.get_amplitudes(other._cpp))
        if isinstance(other, StateAtom):
            return self._cpp.get_amplitudes(other._cpp).toarray().flatten()
        if isinstance(other, BasisAtom):
            return self._cpp.get_amplitudes(other._cpp)
        raise TypeError(f"Incompatible types: {type(other)=}; {type(self)=}")

    @overload
    def get_overlaps(self, other: KetAtom | StateAtom) -> NDArray: ...

    @overload
    def get_overlaps(self, other: Self) -> csr_matrix: ...

    def get_overlaps(self, other: KetAtom | StateAtom | Self) -> NDArray | csr_matrix:
        if isinstance(other, KetAtom):
            return np.array(self._cpp.get_overlaps(other._cpp))
        if isinstance(other, StateAtom):
            return self._cpp.get_overlaps(other._cpp).toarray().flatten()
        if isinstance(other, BasisAtom):
            return self._cpp.get_overlaps(other._cpp)
        raise TypeError(f"Incompatible types: {type(other)=}; {type(self)=}")

    @overload
    def get_matrix_elements(
        self, other: KetAtom | StateAtom, operator: OperatorType, q: int, unit: None = None
    ) -> PintArray: ...

    @overload
    def get_matrix_elements(self, other: KetAtom | StateAtom, operator: OperatorType, q: int, unit: str) -> NDArray: ...

    @overload
    def get_matrix_elements(self, other: Self, operator: OperatorType, q: int, unit: None = None) -> PintSparse: ...

    @overload
    def get_matrix_elements(self, other: Self, operator: OperatorType, q: int, unit: str) -> csr_matrix: ...

    def get_matrix_elements(
        self, other: KetAtom | StateAtom | Self, operator: OperatorType, q: int, unit: str | None = None
    ) -> NDArray | PintArray | csr_matrix | PintSparse:
        cpp_op = get_cpp_operator_type(operator)

        matrix_elements_au: NDArray
        if isinstance(other, KetAtom):
            matrix_elements_au = np.array(self._cpp.get_matrix_elements(other._cpp, cpp_op, q))
            return QuantityArray.convert_au_to_user(matrix_elements_au, operator, unit)
        if isinstance(other, StateAtom):
            matrix_elements_au = self._cpp.get_matrix_elements(other._cpp, cpp_op, q).toarray().flatten()
            return QuantityArray.convert_au_to_user(matrix_elements_au, operator, unit)
        if isinstance(other, BasisAtom):
            matrix_elements_sparse_au = self._cpp.get_matrix_elements(other._cpp, cpp_op, q)
            return QuantitySparse.convert_au_to_user(matrix_elements_sparse_au, operator, unit)
        raise TypeError(f"Unknown type: {type(other)=}")


class BasisAtomReal(BasisAtom):
    _cpp: _backend.BasisAtomReal  # type: ignore [assignment]
    _cpp_creator = _backend.BasisAtomCreatorReal  # type: ignore [assignment]
    _ket_class = KetAtom
    _state_class = StateAtomReal
