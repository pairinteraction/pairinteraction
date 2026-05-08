# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any, TypeAlias, TypeGuard, cast, overload

import numpy as np
from typing_extensions import Self

from pairinteraction import _backend
from pairinteraction.basis.basis_atom import BasisAtom
from pairinteraction.basis.basis_base import BasisBase
from pairinteraction.enums import OperatorType, Parity, get_cpp_operator_type, get_cpp_parity
from pairinteraction.ket import KetPair, KetPairReal, is_ket_atom_tuple
from pairinteraction.state import StatePair, StatePairReal
from pairinteraction.state.state_pair import is_state_atom_tuple
from pairinteraction.units import QuantityArray, QuantityScalar, QuantitySparse

if TYPE_CHECKING:
    from collections.abc import Sequence

    from scipy.sparse import csr_matrix

    from pairinteraction.basis.basis_base import UnionCPPBasis
    from pairinteraction.enums import OperatorType, Parity
    from pairinteraction.ket import KetAtomTuple, KetPairLike
    from pairinteraction.state.state_pair import StatePairLike
    from pairinteraction.system import SystemAtom
    from pairinteraction.units import NDArray, PintArray, PintFloat, PintSparse

BasisPairLike: TypeAlias = "BasisPair | tuple[BasisAtom, BasisAtom] | Sequence[BasisAtom]"

logger = logging.getLogger(__name__)


def is_basis_pair_like(obj: Any) -> TypeGuard[BasisPairLike]:
    return isinstance(obj, BasisPair) or is_basis_atom_tuple(obj)


def is_basis_atom_tuple(obj: Any) -> TypeGuard[tuple[BasisAtom, BasisAtom]]:
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
        >>> import pairinteraction as pi
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
        BasisPair(|Rb:59,S_1/2,-1/2; Rb:61,S_1/2,-1/2⟩ ... |Rb:58,F_7/2,7/2; Rb:59,S_1/2,1/2⟩)

    """

    _cpp: _backend.BasisPairComplex
    _cpp_creator = _backend.BasisPairCreatorComplex
    _ket_class = KetPair
    _state_class = StatePair

    system_atoms: tuple[SystemAtom, SystemAtom]
    """The two SystemAtom objects, from which the BasisPair is build."""

    def __init__(
        self,
        systems: Sequence[SystemAtom],
        m: tuple[float, float] | None = None,
        product_of_parities: Parity | None = None,
        energy: tuple[float, float] | tuple[PintFloat, PintFloat] | None = None,
        energy_unit: str | None = None,
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

        self.system_atoms = tuple(systems)  # type: ignore [assignment]

        self._post_init()

    @classmethod
    def from_ket_atoms(
        cls: type[Self],
        ket_atom_tuples: KetAtomTuple | Sequence[KetAtomTuple],
        system_atoms: Sequence[SystemAtom],
        delta_m: float | None = None,
        product_of_parities: Parity | None = None,
        delta_energy: float | PintFloat | None = None,
        delta_energy_unit: str | None = None,
        number_of_kets: int | None = None,
    ) -> Self:
        """Create a BasisPair from one or more pairs of kets with optional energy/m windows.

        Each element of ``ket_atom_tuples`` is a pair ``(ket1, ket2)`` of
        :class:`~pairinteraction.KetAtom` objects.

        Currently a single big basis including all kets for the quantum numbers
        from min_value - delta to max_value + delta is returned.
        In the future this might change to return a basis including all states around the given kets +/- delta,
        but not necessarily all states between the given kets.

        You can either give a ``delta_energy`` to restrict the basis size in energy,
        or give an (approximate) number_of_kets, which is used to construct a basis centered around
        the provided ket pairs with approximately that many kets. If there are multiple states with the same energy,
        the actual number of kets may be a bit higher than the specified number.

        Args:
            ket_atom_tuples: A single pair ``(ket1, ket2)`` of
                :class:`~pairinteraction.KetAtom` objects, or a list of such pairs.
                Must not be empty.
            system_atoms: A collection of exactly two diagonalized
                :class:`~pairinteraction.SystemAtom` objects, one per atom.
            delta_m: Half-width of the total magnetic quantum number window
                ``m = m1 + m2``. Default None means no m restriction.
            product_of_parities: Restrict to pair states with this product of
                single-atom parities. Default None means no parity restriction.
            delta_energy: Half-width of the energy window. Mutually exclusive with
                ``number_of_kets``. Default None means no energy restriction.
            delta_energy_unit: Unit for ``delta_energy`` and the pair energies
                (e.g. ``"GHz"``). Default None means pint quantities are used.
            number_of_kets: Target number of pair kets to include. The method
                keeps the ``number_of_kets`` states closest in energy to the
                reference ket pairs. Mutually exclusive with ``delta_energy``.
                Default None means no count restriction.

        Returns:
            A new :class:`BasisPair` whose energy (and optionally m) range is
            determined by the provided ket pairs and the chosen restrictions.

        Examples:
            >>> import pairinteraction as pi
            >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
            >>> basis_atom = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2))
            >>> system = pi.SystemAtom(basis_atom).diagonalize()
            >>> pair_basis = pi.BasisPair.from_ket_atoms(
            ...     (ket, ket), [system, system], delta_energy=3, delta_energy_unit="GHz"
            ... )

        """
        if number_of_kets is not None:
            if delta_energy is not None:
                raise ValueError("Cannot specify both number_of_kets and delta_energy. Please choose one of them.")
            if number_of_kets <= 0:
                raise ValueError("number_of_kets must be positive.")
            logger.warning(
                "Use number_of_kets with caution and only if you know what you are doing! "
                "It might lead to unexpected results. "
                "E.g. loosening other restrictions while keeping number_of_kets fixed can lead to worse results!"
            )

        if len(ket_atom_tuples) == 0:
            raise ValueError("ket_atom_tuples must not be empty.")
        if is_ket_atom_tuple(ket_atom_tuples):
            ket_atom_tuples = [ket_atom_tuples]
        if not all(is_ket_atom_tuple(t) for t in ket_atom_tuples):
            raise ValueError("ket_atom_tuples must be a KetAtomTuple or a list of KetAtomTuple.")
        ket_atom_tuples = cast("Sequence[KetAtomTuple]", ket_atom_tuples)

        energy_range = None
        if delta_energy is not None:
            pair_energies = [
                sum(
                    system.get_corresponding_energy(ket, delta_energy_unit)
                    for ket, system in zip(ket_atoms, system_atoms, strict=True)
                )
                for ket_atoms in ket_atom_tuples
            ]
            energy_range = (min(pair_energies) - delta_energy, max(pair_energies) + delta_energy)

        m_range = None
        if delta_m is not None:
            pair_ms = [sum(ket.m for ket in ket_atoms) for ket_atoms in ket_atom_tuples]
            m_range = (min(pair_ms) - delta_m, max(pair_ms) + delta_m)

        basis_pair = cls(
            system_atoms,
            m=m_range,
            product_of_parities=product_of_parities,
            energy=energy_range,
            energy_unit=delta_energy_unit,
        )
        if number_of_kets is None or number_of_kets >= basis_pair.number_of_kets:
            return basis_pair

        pair_energies_au = [
            sum(
                system.get_corresponding_energy(ket, "hartree")
                for ket, system in zip(ket_atoms, system_atoms, strict=True)
            )
            for ket_atoms in ket_atom_tuples
        ]
        min_energy_au, max_energy_au = min(pair_energies_au), max(pair_energies_au)

        cpp_kets = basis_pair._cpp.get_kets()
        all_energies_au = np.array([ket_cpp.get_energy() for ket_cpp in cpp_kets])
        deltas = np.maximum(np.maximum(min_energy_au - all_energies_au, all_energies_au - max_energy_au), 0)
        delta_energy = float(np.partition(deltas, number_of_kets - 1)[number_of_kets - 1]) + 1e-10

        return cls(
            system_atoms,
            m=m_range,
            product_of_parities=product_of_parities,
            energy=(min_energy_au - delta_energy, max_energy_au + delta_energy),
            energy_unit="hartree",
        )

    @classmethod
    def _from_cpp_object(cls: type[Self], cpp_obj: UnionCPPBasis, system_atoms: tuple[SystemAtom, SystemAtom]) -> Self:
        obj = super()._from_cpp_object(cpp_obj)
        obj.system_atoms = system_atoms
        return obj

    def get_corresponding_state(self, ket: KetPairLike) -> StatePair:  # type: ignore [override]
        state_index = self.get_corresponding_state_index(ket)
        return self.get_state(state_index)

    def get_corresponding_state_index(self, ket: KetPairLike) -> int:  # type: ignore [override]
        if isinstance(ket, KetPair):
            return super().get_corresponding_state_index(ket)
        if is_ket_atom_tuple(ket):
            overlaps = self.get_overlaps(ket)
            id_max = np.argmax(overlaps)
            if overlaps[id_max] < 0.51:
                raise ValueError(
                    "The provided ket pair does not correspond well to any state in the basis "
                    f"(max overlap={overlaps[id_max]:.3f})."
                )
            return int(id_max)
        raise TypeError(f"Unknown type: {type(ket)=}")

    @overload
    def get_amplitudes(self, other: KetPairLike | StatePairLike) -> NDArray: ...

    @overload
    def get_amplitudes(self, other: BasisPairLike) -> csr_matrix: ...

    def get_amplitudes(self, other: KetPairLike | StatePairLike | BasisPairLike) -> NDArray | csr_matrix:
        # KetPair like
        if isinstance(other, KetPair):
            return np.array(self._cpp.get_amplitudes(other._cpp))
        if is_ket_atom_tuple(other):
            ket_cpp = (other[0]._cpp, other[1]._cpp)
            return np.array(self._cpp.get_amplitudes(*ket_cpp))

        # StatePair
        if isinstance(other, StatePair):
            return self._cpp.get_amplitudes(other._cpp).toarray().flatten()
        if is_state_atom_tuple(other):
            state_cpp = (other[0]._cpp, other[1]._cpp)
            return self._cpp.get_amplitudes(*state_cpp).toarray().flatten()

        # BasisPair like
        if isinstance(other, BasisPair):
            return self._cpp.get_amplitudes(other._cpp)
        if is_basis_atom_tuple(other):
            basis_cpp = (other[0]._cpp, other[1]._cpp)
            return self._cpp.get_amplitudes(*basis_cpp)

        raise TypeError(f"Unknown type: {type(other)=}")

    @overload
    def get_overlaps(self, other: KetPairLike | StatePairLike) -> NDArray: ...

    @overload
    def get_overlaps(self, other: BasisPairLike) -> csr_matrix: ...

    def get_overlaps(self, other: KetPairLike | StatePairLike | BasisPairLike) -> NDArray | csr_matrix:
        # KetPair like
        if isinstance(other, KetPair):
            return np.array(self._cpp.get_overlaps(other._cpp))
        if is_ket_atom_tuple(other):
            ket_cpp = (other[0]._cpp, other[1]._cpp)
            return np.array(self._cpp.get_overlaps(*ket_cpp))

        # StatePair
        if isinstance(other, StatePair):
            return self._cpp.get_overlaps(other._cpp).toarray().flatten()
        if is_state_atom_tuple(other):
            state_cpp = (other[0]._cpp, other[1]._cpp)
            return self._cpp.get_overlaps(*state_cpp).toarray().flatten()

        # BasisPair like
        if isinstance(other, BasisPair):
            return self._cpp.get_overlaps(other._cpp)
        if is_basis_atom_tuple(other):
            basis_cpp = (other[0]._cpp, other[1]._cpp)
            return self._cpp.get_overlaps(*basis_cpp)

        raise TypeError(f"Unknown type: {type(other)=}")

    @overload
    def get_matrix_elements(
        self,
        other: KetPairLike | StatePairLike,
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        unit: None = None,
    ) -> PintArray: ...

    @overload
    def get_matrix_elements(
        self,
        other: KetPairLike | StatePairLike,
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        unit: str,
    ) -> NDArray: ...

    @overload
    def get_matrix_elements(
        self,
        other: BasisPairLike,
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        unit: None = None,
    ) -> PintSparse: ...

    @overload
    def get_matrix_elements(
        self,
        other: BasisPairLike,
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        unit: str,
    ) -> csr_matrix: ...

    def get_matrix_elements(
        self,
        other: KetPairLike | StatePairLike | BasisPairLike,
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        unit: str | None = None,
    ) -> NDArray | PintArray | csr_matrix | PintSparse:
        operators_cpp = (get_cpp_operator_type(operators[0]), get_cpp_operator_type(operators[1]))
        matrix_elements_au: NDArray

        # KetPair like
        if isinstance(other, KetPair):
            matrix_elements_au = np.array(self._cpp.get_matrix_elements(other._cpp, *operators_cpp, *qs))
            return QuantityArray.convert_au_to_user(matrix_elements_au, operators, unit)
        if is_ket_atom_tuple(other):
            ket_cpp = (other[0]._cpp, other[1]._cpp)
            matrix_elements_au = np.array(self._cpp.get_matrix_elements(*ket_cpp, *operators_cpp, *qs))
            return QuantityArray.convert_au_to_user(matrix_elements_au, operators, unit)

        # StatePair like
        if isinstance(other, StatePair):
            matrix_elements_au = self._cpp.get_matrix_elements(other._cpp, *operators_cpp, *qs).toarray().flatten()
            return QuantityArray.convert_au_to_user(matrix_elements_au, operators, unit)
        if is_state_atom_tuple(other):
            state_cpp = (other[0]._cpp, other[1]._cpp)
            matrix_elements_au = self._cpp.get_matrix_elements(*state_cpp, *operators_cpp, *qs).toarray().flatten()
            return QuantityArray.convert_au_to_user(matrix_elements_au, operators, unit)

        # BasisPair like
        if isinstance(other, BasisPair):
            matrix_elements_sparse_au = self._cpp.get_matrix_elements(other._cpp, *operators_cpp, *qs)
            return QuantitySparse.convert_au_to_user(matrix_elements_sparse_au, operators, unit)
        if is_basis_atom_tuple(other):
            basis_cpp = (other[0]._cpp, other[1]._cpp)
            matrix_elements_sparse_au = self._cpp.get_matrix_elements(*basis_cpp, *operators_cpp, *qs)
            return QuantitySparse.convert_au_to_user(matrix_elements_sparse_au, operators, unit)

        raise TypeError(f"Unknown type: {type(other)=}")


class BasisPairReal(BasisPair):
    _cpp: _backend.BasisPairReal  # type: ignore [assignment]
    _cpp_creator = _backend.BasisPairCreatorReal  # type: ignore [assignment]
    _ket_class = KetPairReal
    _state_class = StatePairReal
