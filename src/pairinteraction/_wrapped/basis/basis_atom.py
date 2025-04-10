# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Any, ClassVar, Optional, TypeVar, Union, overload

import numpy as np

from pairinteraction import _backend
from pairinteraction._wrapped.basis.basis import BasisBase
from pairinteraction._wrapped.database.database import Database
from pairinteraction._wrapped.enums import OperatorType, Parity, get_cpp_operator_type, get_cpp_parity
from pairinteraction._wrapped.ket.ket_atom import KetAtom
from pairinteraction._wrapped.state.state_atom import StateAtom, StateAtomComplex, StateAtomReal
from pairinteraction.units import QuantityArray, QuantityScalar, QuantitySparse

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix
    from typing_extensions import Self

    from pairinteraction.units import NDArray, PintArray, PintFloat, PintSparse

StateType = TypeVar("StateType", bound="StateAtom[Any]", covariant=True)
UnionCPPBasisAtom = Union[_backend.BasisAtomReal, _backend.BasisAtomComplex]
UnionTypeCPPBasisAtomCreator = Union[type[_backend.BasisAtomCreatorReal], type[_backend.BasisAtomCreatorComplex]]


class BasisAtom(BasisBase[KetAtom, StateType]):
    """Basis for a single atom.

    Add all KetAtom objects that match the given quantum numbers to the basis.
    The initial coefficients matrix is a unit matrix, i.e. the first basis state is the first ket, etc.
    The BasisAtom coefficients matrix will always be square,
    i.e. the number of kets is equal to the number of states.

    Examples:
        >>> import pairinteraction.real as pi
        >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
        >>> energy_min, energy_max = ket.get_energy(unit="GHz") - 100, ket.get_energy(unit="GHz") + 100
        >>> basis = pi.BasisAtom("Rb", n=(57, 63), l=(0, 3), energy=(energy_min, energy_max), energy_unit="GHz")
        >>> print(basis)
        BasisAtom(n=(57, 63), l=(0, 3), energy=(1008911.9216, 1009111.9216), energy_unit=GHz)

    """

    _cpp: UnionCPPBasisAtom
    _cpp_creator: ClassVar[UnionTypeCPPBasisAtomCreator]

    def __init__(  # noqa: C901, PLR0912
        self,
        species: str,
        n: Optional[tuple[int, int]] = None,
        nu: Optional[tuple[float, float]] = None,
        l: Optional[tuple[float, float]] = None,
        s: Optional[tuple[float, float]] = None,
        j: Optional[tuple[float, float]] = None,
        f: Optional[tuple[float, float]] = None,
        m: Optional[tuple[float, float]] = None,
        energy: Union[tuple[float, float], tuple["PintFloat", "PintFloat"], None] = None,
        energy_unit: Optional[str] = None,
        parity: Optional[Parity] = None,
        database: Optional[Database] = None,
        additional_kets: Optional[list[KetAtom]] = None,
    ) -> None:
        """Create a basis for a single atom.

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
                Default None, i.e. energy is provided as pint object.
            parity: The parity of the states to consider. Default None, i.e. add all available states.
            database: Which database to use. Default None, i.e. use the global database instance.
            additional_kets: List of additional kets to add to the basis. Default None.

        """
        creator = self._cpp_creator()
        creator.set_species(species)
        self.species = species

        self._qns: dict[str, tuple[float, float]] = {}
        if n is not None:
            if not all(isinstance(x, int) or x.is_integer() for x in n):
                raise ValueError("Quantum numbers n must be integers.")
            n = (int(n[0]), int(n[1]))
            creator.restrict_quantum_number_n(*n)
            self._qns["n"] = n
        if nu is not None:
            creator.restrict_quantum_number_nu(*nu)
            self._qns["nu"] = nu
        if s is not None:
            creator.restrict_quantum_number_s(*s)
            self._qns["s"] = s
        if l is not None:
            creator.restrict_quantum_number_l(*l)
            self._qns["l"] = l
        if j is not None:
            self._qns["j"] = j
            creator.restrict_quantum_number_j(*j)
        if f is not None:
            creator.restrict_quantum_number_f(*f)
            self._qns["f"] = f
        if m is not None:
            creator.restrict_quantum_number_m(*m)
            self._qns["m"] = m
        if parity is not None:
            creator.restrict_parity(get_cpp_parity(parity))
        if energy is not None:
            min_energy_au = QuantityScalar.from_pint_or_unit(energy[0], energy_unit, "energy").to_base_unit()
            max_energy_au = QuantityScalar.from_pint_or_unit(energy[1], energy_unit, "energy").to_base_unit()
            creator.restrict_energy(min_energy_au, max_energy_au)
            self._energy = energy
            self._energy_unit = energy_unit
        if database is None:
            if Database.get_global_database() is None:
                Database.initialize_global_database()
            database = Database.get_global_database()
        if additional_kets is not None:
            for ket in additional_kets:
                creator.append_ket(ket._cpp)
            self._additional_kets = additional_kets
        self._database = database
        self._cpp = creator.create(database._cpp)

    def __repr__(self) -> str:
        args = ""
        if hasattr(self, "_qns"):
            args += ", ".join(f"{k}={v}" for k, v in self._qns.items())
        if hasattr(self, "_energy"):
            args += f", energy=({self._energy[0]:.4f}, {self._energy[1]:.4f}), energy_unit={self._energy_unit}"
        if hasattr(self, "_additional_kets"):
            args += f", additional_kets={self._additional_kets}"
        if len(args) == 0:
            return super().__repr__()
        return f"{type(self).__name__}({args})"

    @property
    def database(self) -> Database:
        """The database used for this object."""
        return self._database

    @overload
    def get_amplitudes(self, other: Union[KetAtom, StateAtom[Any]]) -> "NDArray": ...

    @overload
    def get_amplitudes(self, other: "Self") -> "csr_matrix": ...

    def get_amplitudes(self, other: Union[KetAtom, StateAtom[Any], "Self"]) -> Union["NDArray", "csr_matrix"]:
        if isinstance(other, KetAtom):
            return np.array(self._cpp.get_amplitudes(other._cpp))
        if isinstance(self, BasisAtomReal) and isinstance(other, StateAtomReal):
            return self._cpp.get_amplitudes(other._basis._cpp).toarray().flatten()
        if isinstance(self, BasisAtomComplex) and isinstance(other, StateAtomComplex):
            return self._cpp.get_amplitudes(other._basis._cpp).toarray().flatten()
        if isinstance(self, BasisAtomReal) and isinstance(other, BasisAtomReal):
            return self._cpp.get_amplitudes(other._cpp)
        if isinstance(self, BasisAtomComplex) and isinstance(other, BasisAtomComplex):
            return self._cpp.get_amplitudes(other._cpp)
        raise TypeError(f"Incompatible types: {type(other)=}; {type(self)=}")

    @overload
    def get_overlaps(self, other: Union[KetAtom, StateAtom[Any]]) -> "NDArray": ...

    @overload
    def get_overlaps(self, other: "Self") -> "csr_matrix": ...

    def get_overlaps(self, other: Union[KetAtom, StateAtom[Any], "Self"]) -> Union["NDArray", "csr_matrix"]:
        if isinstance(other, KetAtom):
            return np.array(self._cpp.get_overlaps(other._cpp))
        if isinstance(self, BasisAtomReal) and isinstance(other, StateAtomReal):
            return self._cpp.get_overlaps(other._basis._cpp).toarray().flatten()
        if isinstance(self, BasisAtomComplex) and isinstance(other, StateAtomComplex):
            return self._cpp.get_overlaps(other._basis._cpp).toarray().flatten()
        if isinstance(self, BasisAtomReal) and isinstance(other, BasisAtomReal):
            return self._cpp.get_overlaps(other._cpp)
        if isinstance(self, BasisAtomComplex) and isinstance(other, BasisAtomComplex):
            return self._cpp.get_overlaps(other._cpp)
        raise TypeError(f"Incompatible types: {type(other)=}; {type(self)=}")

    @overload
    def get_matrix_elements(
        self, other: Union[KetAtom, StateAtom[Any]], operator: OperatorType, q: int, unit: None = None
    ) -> "PintArray": ...

    @overload
    def get_matrix_elements(
        self, other: Union[KetAtom, StateAtom[Any]], operator: OperatorType, q: int, unit: str
    ) -> "NDArray": ...

    @overload
    def get_matrix_elements(self, other: "Self", operator: OperatorType, q: int, unit: None = None) -> "PintSparse": ...  # type: ignore [type-var] # see PintSparse

    @overload
    def get_matrix_elements(self, other: "Self", operator: OperatorType, q: int, unit: str) -> "csr_matrix": ...

    def get_matrix_elements(
        self, other: Union[KetAtom, StateAtom[Any], "Self"], operator: OperatorType, q: int, unit: Optional[str] = None
    ) -> Union["NDArray", "PintArray", "csr_matrix", "PintSparse"]:
        cpp_op = get_cpp_operator_type(operator)

        matrix_elements_au: Union[NDArray, csr_matrix]
        matrix_elements: Union[QuantityArray, QuantitySparse]
        if isinstance(other, KetAtom):
            matrix_elements_au = np.array(self._cpp.get_matrix_elements(other._cpp, cpp_op, q))
            matrix_elements = QuantityArray.from_base_unit(matrix_elements_au, operator)
        elif isinstance(other, StateAtom):
            if isinstance(self, BasisAtomReal) and isinstance(other, StateAtomReal):
                matrix_elements_au = self._cpp.get_matrix_elements(other._basis._cpp, cpp_op, q).toarray().flatten()
            elif isinstance(self, BasisAtomComplex) and isinstance(other, StateAtomComplex):
                matrix_elements_au = self._cpp.get_matrix_elements(other._basis._cpp, cpp_op, q).toarray().flatten()
            else:
                raise TypeError(f"Incompatible types: {type(other)=}; {type(self)=}")
            matrix_elements = QuantityArray.from_base_unit(matrix_elements_au, operator)
        elif isinstance(other, BasisAtom):
            if isinstance(self, BasisAtomReal) and isinstance(other, BasisAtomReal):
                matrix_elements_au = self._cpp.get_matrix_elements(other._cpp, cpp_op, q)
            elif isinstance(self, BasisAtomComplex) and isinstance(other, BasisAtomComplex):
                matrix_elements_au = self._cpp.get_matrix_elements(other._cpp, cpp_op, q)
            else:
                raise TypeError(f"Incompatible types: {type(other)=}; {type(self)=}")
            matrix_elements = QuantitySparse.from_base_unit(matrix_elements_au, operator)
        else:
            raise TypeError(f"Unknown type: {type(other)=}")

        return matrix_elements.to_pint_or_unit(unit)


class BasisAtomReal(BasisAtom[StateAtomReal]):
    _cpp: _backend.BasisAtomReal
    _cpp_creator = _backend.BasisAtomCreatorReal
    _TypeKet = KetAtom
    _TypeState = StateAtomReal


class BasisAtomComplex(BasisAtom[StateAtomComplex]):
    _cpp: _backend.BasisAtomComplex
    _cpp_creator = _backend.BasisAtomCreatorComplex
    _TypeKet = KetAtom
    _TypeState = StateAtomComplex
