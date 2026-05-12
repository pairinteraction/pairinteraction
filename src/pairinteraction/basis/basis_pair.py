# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any, TypeAlias, TypeGuard, cast, overload

import numpy as np
from typing_extensions import Self, deprecated

from pairinteraction import _backend
from pairinteraction.basis.basis_atom import BasisAtom
from pairinteraction.basis.basis_base import BasisBase
from pairinteraction.enums import OperatorType, Parity, get_cpp_operator_type, get_cpp_parity
from pairinteraction.ket import KetPair, KetPairReal, is_ket_atom_tuple
from pairinteraction.ket.ket_pair import get_ketpairlike_energy, get_ketpairlike_m, is_ket_pair_like
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
    and which parities under inversion and permutation should be used for symmetrization.
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
        system_atoms: Sequence[SystemAtom],
        m: tuple[float, float] | None = None,
        parity_under_inversion: Parity | None = None,
        parity_under_permutation: Parity | None = None,
        energy: tuple[float, float] | tuple[PintFloat, PintFloat] | None = None,
        energy_unit: str | None = None,
    ) -> None:
        """Create a basis for a pair of atoms.

        Args:
            system_atoms: tuple of two SystemAtom objects, which define the two atoms, from which the BasisPair is build
                Both system_atoms have to be diagonalized before creating the BasisPair.
            m: tuple of (min, max) values for the total magnetic quantum number m of the pair state.
                Default None, i.e. no restriction.
            parity_under_inversion: Restrict to pair states with this parity under inversion.
                Default None, i.e. do not apply inversion symmetrization.
            parity_under_permutation: Restrict to pair states with this parity under permutation.
                Default None, i.e. do not apply permutation symmetrization.
            energy: tuple of (min, max) value for the pair energy. Default None, i.e. add all available states.
            energy_unit: In which unit the energy values are given, e.g. "GHz".
                Default None, i.e. energy is provided as pint object.

        """
        assert len(system_atoms) == 2, "BasisPair requires exactly two SystemAtom objects."
        creator = self._cpp_creator()
        for system in system_atoms:
            creator.add(system._cpp)
        if m is not None:
            creator.restrict_quantum_number_m(*m)
        if parity_under_inversion is not None:
            creator.restrict_parity_under_inversion(get_cpp_parity(parity_under_inversion))
        if parity_under_permutation is not None:
            creator.restrict_parity_under_permutation(get_cpp_parity(parity_under_permutation))
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

        self.system_atoms = tuple(system_atoms)  # type: ignore [assignment]

        self._post_init()

    @classmethod
    def from_kets(  # noqa: C901
        cls: type[Self],
        kets: KetPairLike | Sequence[KetPairLike],
        system_atoms: Sequence[SystemAtom],
        delta_m: float | None = None,
        parity_under_inversion: Parity | None = None,
        parity_under_permutation: Parity | None = None,
        delta_energy: float | PintFloat | None = None,
        delta_energy_unit: str | None = None,
        number_of_kets: int | None = None,
        *,
        warn_number_of_kets: bool = True,
    ) -> Self:
        """Create a BasisPair from one or more pairs of kets with optional energy/m windows.

        Currently a single big basis including all kets for the quantum numbers
        from min_value - delta to max_value + delta is returned.
        In the future this might change to return a basis including all states around the given kets +/- delta,
        but not necessarily all states between the given kets.

        You can either give a ``delta_energy`` to restrict the basis size in energy,
        or give an (approximate) number_of_kets, which is used to construct a basis centered around
        the provided ket pairs with approximately that many kets. If there are multiple states with the same energy,
        the actual number of kets may be a bit higher than the specified number.

        Args:
            kets: A single ket pair (given either as tuple ``(ket1, ket2)`` of
                :class:`~pairinteraction.KetAtom` objects or as a :class:`~pairinteraction.KetPair` object)
                or a list of ket pairs.
                Must not be empty.
            system_atoms: A collection of exactly two diagonalized
                :class:`~pairinteraction.SystemAtom` objects, one per atom.
            delta_m: Half-width of the total magnetic quantum number window
                ``m = m1 + m2``. Default None means no m restriction.
            parity_under_inversion: Restrict to pair states with this parity under inversion.
                Default None means no inversion symmetrization.
            parity_under_permutation: Restrict to pair states with this parity under permutation.
                Default None means no permutation symmetrization.
            delta_energy: Half-width of the energy window. Mutually exclusive with
                ``number_of_kets``. Default None means no energy restriction.
            delta_energy_unit: Unit for ``delta_energy`` and the pair energies
                (e.g. ``"GHz"``). Default None means pint quantities are used.
            number_of_kets: Target number of pair kets to include. The method
                keeps the ``number_of_kets`` states closest in energy to the
                reference ket pairs. Mutually exclusive with ``delta_energy``.
                Default None means no count restriction.
            warn_number_of_kets: Don't warn about the possible issues with using number_of_kets.

        Returns:
            A new :class:`BasisPair` whose energy (and optionally m) range is
            determined by the provided ket pairs and the chosen restrictions.

        Examples:
            >>> import pairinteraction as pi
            >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
            >>> basis_atom = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2))
            >>> system = pi.SystemAtom(basis_atom).diagonalize()
            >>> pair_basis = pi.BasisPair.from_kets(
            ...     (ket, ket), [system, system], delta_energy=3, delta_energy_unit="GHz"
            ... )

        """
        assert len(system_atoms) == 2, "BasisPair requires exactly two SystemAtom objects."

        if number_of_kets is not None:
            if delta_energy is not None:
                raise ValueError("Cannot specify both number_of_kets and delta_energy. Please choose one of them.")
            if number_of_kets <= 0:
                raise ValueError("number_of_kets must be positive.")
            if warn_number_of_kets:
                logger.warning(
                    "Use number_of_kets with caution and only if you know what you are doing! "
                    "It might lead to unexpected results. "
                    "E.g. loosening other restrictions while keeping number_of_kets fixed can lead to worse results!"
                )

        if is_ket_atom_tuple(kets) or isinstance(kets, KetPair):
            kets = [kets]
        if not all(is_ket_pair_like(t) for t in kets):
            raise ValueError("kets must be a KetPairLike or a list of KetPairLike.")
        if len(kets) == 0:
            raise ValueError("kets must not be empty.")
        kets = cast("Sequence[KetAtomTuple | KetPair]", kets)

        energy_range = None
        if delta_energy is not None:
            pair_energies = [get_ketpairlike_energy(ket, system_atoms, delta_energy_unit) for ket in kets]
            energy_range = (min(pair_energies) - delta_energy, max(pair_energies) + delta_energy)

        m_range = None
        if delta_m is not None:
            pair_ms = [get_ketpairlike_m(ket) for ket in kets]
            m_range = (min(pair_ms) - delta_m, max(pair_ms) + delta_m)

        basis_pair = cls(
            system_atoms,
            m=m_range,
            parity_under_inversion=parity_under_inversion,
            parity_under_permutation=parity_under_permutation,
            energy=energy_range,
            energy_unit=delta_energy_unit,
        )
        if number_of_kets is None or number_of_kets >= basis_pair.number_of_kets:
            return basis_pair

        pair_energies_au = [get_ketpairlike_energy(ket, system_atoms, "hartree") for ket in kets]
        min_energy_au, max_energy_au = min(pair_energies_au), max(pair_energies_au)

        cpp_kets = basis_pair._cpp.get_kets()
        all_energies_au = np.array([ket_cpp.get_energy() for ket_cpp in cpp_kets])
        deltas = np.maximum(np.maximum(min_energy_au - all_energies_au, all_energies_au - max_energy_au), 0)
        delta_energy = float(np.partition(deltas, number_of_kets - 1)[number_of_kets - 1]) + 1e-10

        return cls(
            system_atoms,
            m=m_range,
            parity_under_inversion=parity_under_inversion,
            parity_under_permutation=parity_under_permutation,
            energy=(min_energy_au - delta_energy, max_energy_au + delta_energy),
            energy_unit="hartree",
        )

    @classmethod
    @deprecated("Use `BasisPair.from_kets` instead.")
    def from_ket_atoms(
        cls: type[Self],
        ket_atom_tuples: KetAtomTuple | Sequence[KetAtomTuple],
        system_atoms: Sequence[SystemAtom],
        delta_m: float | None = None,
        parity_under_inversion: Parity | None = None,
        parity_under_permutation: Parity | None = None,
        delta_energy: float | PintFloat | None = None,
        delta_energy_unit: str | None = None,
        number_of_kets: int | None = None,
    ) -> Self:
        return cls.from_kets(
            ket_atom_tuples,
            system_atoms,
            delta_m=delta_m,
            parity_under_inversion=parity_under_inversion,
            parity_under_permutation=parity_under_permutation,
            delta_energy=delta_energy,
            delta_energy_unit=delta_energy_unit,
            number_of_kets=number_of_kets,
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
            return self._cpp.get_amplitudes(other._cpp).toarray().ravel()
        if is_state_atom_tuple(other):
            state_cpp = (other[0]._cpp, other[1]._cpp)
            return self._cpp.get_amplitudes(*state_cpp).toarray().ravel()

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
            return self._cpp.get_overlaps(other._cpp).toarray().ravel()
        if is_state_atom_tuple(other):
            state_cpp = (other[0]._cpp, other[1]._cpp)
            return self._cpp.get_overlaps(*state_cpp).toarray().ravel()

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
            matrix_elements_au = self._cpp.get_matrix_elements(other._cpp, *operators_cpp, *qs).toarray().ravel()
            return QuantityArray.convert_au_to_user(matrix_elements_au, operators, unit)
        if is_state_atom_tuple(other):
            state_cpp = (other[0]._cpp, other[1]._cpp)
            matrix_elements_au = self._cpp.get_matrix_elements(*state_cpp, *operators_cpp, *qs).toarray().ravel()
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
