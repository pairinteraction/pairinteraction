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
        additional_kets: list[KetAtom] | None = None,
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
        return self.kets[0].database

    @property
    def species(self) -> str:
        """The atomic species."""
        return self.kets[0].species

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
    def get_matrix_elements(self, other: Self, operator: OperatorType, q: int, unit: None = None) -> PintSparse: ...  # type: ignore [type-var] # see PintSparse

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
