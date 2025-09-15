# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from collections.abc import Collection, Sequence
from typing import TYPE_CHECKING, Any, Optional, Union, overload

import numpy as np
from typing_extensions import TypeGuard

from pairinteraction import _backend
from pairinteraction.basis.basis import BasisBase
from pairinteraction.basis.basis_atom import BasisAtom
from pairinteraction.enums import OperatorType, Parity, get_cpp_operator_type, get_cpp_parity
from pairinteraction.ket import KetPair, KetPairReal, is_ket_atom_tuple
from pairinteraction.state import StatePair, StatePairReal
from pairinteraction.units import QuantityArray, QuantityScalar, QuantitySparse

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix

    from pairinteraction.ket import (
        KetAtom,  # noqa: F401  # required for sphinx for KetPairLike
        KetPairLike,
    )
    from pairinteraction.system import SystemAtom
    from pairinteraction.units import NDArray, PintArray, PintFloat, PintSparse

BasisPairLike = Union["BasisPair", tuple["BasisAtom", "BasisAtom"], Sequence["BasisAtom"]]


def is_basis_pair_like(obj: Any) -> TypeGuard[BasisPairLike]:
    return isinstance(obj, BasisPair) or is_basis_atom_tuple(obj)


def is_basis_atom_tuple(obj: Any) -> TypeGuard[tuple["BasisAtom", "BasisAtom"]]:
    return hasattr(obj, "__len__") and len(obj) == 2 and all(isinstance(x, BasisAtom) for x in obj)


class BasisPair(BasisBase[KetPair, StatePair]):
    """Basis for a pair of atoms.

    Add all product states of the eigenstates of two given SystemAtom objects to the basis,
    which pair energy is within the given energy range.
    You can also specify which total magnetic quantum number m the pair should have (if it is conserved)
    and the product of the parities of the two atoms.
    Due to the possible restrictions of the basis states, the BasisPair coefficients matrix will in general
    not be square but (n x d),
    where n is the number of all involved kets (typically basis1.number_of_kets * basis2.number_of_kets)
    and d is the number of basis states (after applying the restrictions).

    Examples:
        >>> import pairinteraction.real as pi
        >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
        >>> basis = pi.BasisAtom("Rb", n=(58, 63), l=(0, 3))
        >>> system = pi.SystemAtom(basis).set_magnetic_field([0, 0, 1], unit="G").diagonalize()
        >>> pair_energy = 2 * system.get_corresponding_energy(ket, unit="GHz")
        >>> pair_basis = pi.BasisPair(
        ...     [system, system],
        ...     energy=(pair_energy - 3, pair_energy + 3),
        ...     energy_unit="GHz",
        ... )
        >>> print(pair_basis)
        BasisPair(|Rb:59,S_1/2,-1/2; Rb:61,S_1/2,-1/2⟩ ... |~Rb:58,F_7/2,7/2; Rb:59,S_1/2,1/2⟩)

    """

    _cpp: _backend.BasisPairComplex
    _cpp_creator = _backend.BasisPairCreatorComplex
    _ket_class = KetPair
    _state_class = StatePair

    def __init__(
        self,
        systems: Collection["SystemAtom"],
        m: Optional[tuple[float, float]] = None,
        product_of_parities: Optional[Parity] = None,
        energy: Union[tuple[float, float], tuple["PintFloat", "PintFloat"], None] = None,
        energy_unit: Optional[str] = None,
    ) -> None:
        """Create a basis for a pair of atoms.

        Args:
            systems: tuple of two SystemAtom objects, which define the two atoms, from which the BasisPair is build.
                Both systems have to be diagonalized before creating the BasisPair.
            m: tuple of (min, max) values for the total magnetic quantum number m of the pair state.
                Default None, i.e. no restriction.
            product_of_parities: The product parity of the states to consider.
                Default None, i.e. add all available states.
            energy: tuple of (min, max) value for the pair energy. Default None, i.e. add all available states.
            energy_unit: In which unit the energy values are given, e.g. "GHz".
                Default None, i.e. energy is provided as pint object.

        """
        assert len(systems) == 2, "BasisPair requires exactly two SystemAtom objects."
        creator = self._cpp_creator()
        self.system_atoms: tuple[SystemAtom, SystemAtom] = tuple(systems)  # type: ignore [assignment]
        for system in systems:
            creator.add(system._cpp)
        if m is not None:
            creator.restrict_quantum_number_m(*m)
        if product_of_parities is not None:
            creator.restrict_product_of_parities(get_cpp_parity(product_of_parities))
        if energy is not None:
            min_energy_au = QuantityScalar.convert_user_to_au(energy[0], energy_unit, "energy")
            max_energy_au = QuantityScalar.convert_user_to_au(energy[1], energy_unit, "energy")
            # in atomic units all energies should be on the order of -0.5 * Z^2 to 0
            # so choosing some very large values for the limits should be fine
            # (we cant use np.inf here, since this is passed to cpp code)
            min_energy_au = np.clip(min_energy_au, -1e10, 1e10)  # FIXME
            max_energy_au = np.clip(max_energy_au, -1e10, 1e10)  # FIXME
            creator.restrict_energy(min_energy_au, max_energy_au)
        self._cpp = creator.create()

    @overload
    def get_amplitudes(self, other: "KetPairLike") -> "NDArray": ...

    @overload
    def get_amplitudes(self, other: BasisPairLike) -> "csr_matrix": ...

    def get_amplitudes(self, other: Union["KetPairLike", BasisPairLike]) -> Union["NDArray", "csr_matrix"]:
        # KetPair like
        if isinstance(other, KetPair):
            return np.array(self._cpp.get_amplitudes(other._cpp))  # type: ignore [call-overload] # _backend.pyi is incorrect
        if is_ket_atom_tuple(other):
            ket_cpp = (other[0]._cpp, other[1]._cpp)
            return np.array(self._cpp.get_amplitudes(*ket_cpp))

        # BasisPair like
        if isinstance(other, BasisPair):
            return self._cpp.get_amplitudes(other._cpp)  # type: ignore [call-overload,no-any-return] # _backend.pyi is incorrect
        if is_basis_atom_tuple(other):
            basis_cpp = (other[0]._cpp, other[1]._cpp)
            return self._cpp.get_amplitudes(*basis_cpp)

        raise TypeError(f"Unknown type: {type(other)=}")

    @overload
    def get_overlaps(self, other: "KetPairLike") -> "NDArray": ...

    @overload
    def get_overlaps(self, other: BasisPairLike) -> "csr_matrix": ...

    def get_overlaps(self, other: Union["KetPairLike", BasisPairLike]) -> Union["NDArray", "csr_matrix"]:
        # KetPair like
        if isinstance(other, KetPair):
            return np.array(self._cpp.get_overlaps(other._cpp))  # type: ignore [call-overload] # _backend.pyi is incorrect
        if is_ket_atom_tuple(other):
            ket_cpp = (other[0]._cpp, other[1]._cpp)
            return np.array(self._cpp.get_overlaps(*ket_cpp))

        # BasisPair like
        if isinstance(other, BasisPair):
            return self._cpp.get_overlaps(other._cpp)  # type: ignore [call-overload,no-any-return] # _backend.pyi is incorrect
        if is_basis_atom_tuple(other):
            basis_cpp = (other[0]._cpp, other[1]._cpp)
            return self._cpp.get_overlaps(*basis_cpp)

        raise TypeError(f"Unknown type: {type(other)=}")

    @overload
    def get_matrix_elements(
        self,
        other: "KetPairLike",
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        unit: None = None,
    ) -> "PintArray": ...

    @overload
    def get_matrix_elements(
        self, other: "KetPairLike", operators: tuple[OperatorType, OperatorType], qs: tuple[int, int], unit: str
    ) -> "NDArray": ...

    @overload
    def get_matrix_elements(
        self,
        other: BasisPairLike,
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        unit: None = None,
    ) -> "PintSparse": ...  # type: ignore [type-var] # see PintSparse

    @overload
    def get_matrix_elements(
        self,
        other: BasisPairLike,
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        unit: str,
    ) -> "csr_matrix": ...

    def get_matrix_elements(
        self,
        other: Union["KetPairLike", BasisPairLike],
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        unit: Optional[str] = None,
    ) -> Union["NDArray", "PintArray", "csr_matrix", "PintSparse"]:
        operators_cpp = (get_cpp_operator_type(operators[0]), get_cpp_operator_type(operators[1]))

        # KetPair like
        if isinstance(other, KetPair):
            matrix_elements_au = np.array(self._cpp.get_matrix_elements(other._cpp, *operators_cpp, *qs))
            return QuantityArray.convert_au_to_user(matrix_elements_au, operators, unit)
        if is_ket_atom_tuple(other):
            ket_cpp = (other[0]._cpp, other[1]._cpp)
            matrix_elements_au = np.array(self._cpp.get_matrix_elements(*ket_cpp, *operators_cpp, *qs))
            return QuantityArray.convert_au_to_user(matrix_elements_au, operators, unit)

        # BasisPair like
        if isinstance(other, BasisPair):
            matrix_elements_au = self._cpp.get_matrix_elements(other._cpp, *operators_cpp, *qs)
            return QuantitySparse.convert_au_to_user(matrix_elements_au, operators, unit)
        if is_basis_atom_tuple(other):
            basis_cpp = (other[0]._cpp, other[1]._cpp)
            matrix_elements_au = self._cpp.get_matrix_elements(*basis_cpp, *operators_cpp, *qs)
            return QuantitySparse.convert_au_to_user(matrix_elements_au, operators, unit)

        raise TypeError(f"Unknown type: {type(other)=}")


class BasisPairReal(BasisPair):
    _cpp: _backend.BasisPairReal
    _cpp_creator = _backend.BasisPairCreatorReal
    _ket_class = KetPairReal
    _state_class = StatePairReal
