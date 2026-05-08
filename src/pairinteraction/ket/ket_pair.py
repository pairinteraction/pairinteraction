# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal, TypeAlias, TypeGuard, overload

from pairinteraction.ket.ket_atom import KetAtom
from pairinteraction.ket.ket_base import KetBase

if TYPE_CHECKING:
    from collections.abc import Sequence

    from pairinteraction import _backend
    from pairinteraction.state import StateAtom
    from pairinteraction.system.system_atom import SystemAtom
    from pairinteraction.units import PintFloat


KetAtomTuple: TypeAlias = "tuple[KetAtom, KetAtom] | Sequence[KetAtom]"
KetPairLike: TypeAlias = "KetPair | KetAtomTuple"


def is_ket_pair_like(obj: Any) -> TypeGuard[KetPairLike]:
    return isinstance(obj, KetPair) or is_ket_atom_tuple(obj)


def is_ket_atom_tuple(obj: Any) -> TypeGuard[KetAtomTuple]:
    return hasattr(obj, "__len__") and len(obj) == 2 and all(isinstance(x, KetAtom) for x in obj)


class KetPair(KetBase):
    """Ket for a pair state of two atoms.

    For pair systems, we choose KetPair object as the product states of the single-atom eigenstates.
    Thus, the Ket pair objects depend on the system and the applied fields.
    Therefore for different pair systems the KetPair objects are not necessarily orthogonal anymore.

    Currently one cannot create a KetPair object directly, but they are used in the background when creating a
    :class:`pairinteraction.BasisPair` object.

    """

    _cpp: _backend.KetPairComplex

    def __init__(self) -> None:
        """Creating a KetPair object directly is not possible."""  # noqa: D401
        raise NotImplementedError("KetPair objects cannot be created directly.")

    def get_label(self, fmt: Literal["raw", "ket", "bra", "detailed"] = "raw", *, max_kets: int = 3) -> str:
        """Label representing the ket pair.

        Args:
            fmt: The format of the label, i.e. whether to return the raw label, or the label in ket or bra notation.
            max_kets: Maximum number of single atom kets to include in the label for each StateAtom.

        Returns:
            A string representation of the ket pair.

        """
        if fmt == "detailed":
            atom_labels = [atom.get_label(max_kets=max_kets) for atom in self.state_atoms]
            return f"({atom_labels[0]}) ⊗ ({atom_labels[1]})"
        return super().get_label(fmt)

    @property
    def state_atoms(self) -> tuple[StateAtom, StateAtom]:
        """Return the state atoms of the ket pair."""
        from pairinteraction.state import StateAtom, StateAtomReal

        _state_atom_class = StateAtomReal if isinstance(self, KetPairReal) else StateAtom

        state_atoms = []
        for atomic_state in self._cpp.get_atomic_states():
            state = _state_atom_class._from_cpp_object(atomic_state)
            state_atoms.append(state)
        return tuple(state_atoms)  # type: ignore [return-value]


class KetPairReal(KetPair):
    _cpp: _backend.KetPairReal  # type: ignore [assignment]


def get_ketpairlike_m(ket: KetPair | KetAtomTuple) -> float:
    if is_ket_atom_tuple(ket):
        m1 = ket[0].m
        m2 = ket[1].m
        return m1 + m2
    if isinstance(ket, KetPair):
        return ket.m
    raise TypeError(f"Unknown type: {type(ket)=}")


@overload
def get_ketpairlike_energy(
    ket: KetPair | KetAtomTuple, system_atoms: Sequence[SystemAtom], unit: None
) -> PintFloat: ...


@overload
def get_ketpairlike_energy(ket: KetPair | KetAtomTuple, system_atoms: Sequence[SystemAtom], unit: str) -> float: ...


def get_ketpairlike_energy(
    ket: KetPair | KetAtomTuple, system_atoms: Sequence[SystemAtom], unit: str | None
) -> float | PintFloat:
    if is_ket_atom_tuple(ket):
        energy1 = system_atoms[0].get_corresponding_energy(ket[0], unit)
        energy2 = system_atoms[1].get_corresponding_energy(ket[1], unit)
        return energy1 + energy2
    if isinstance(ket, KetPair):
        return ket.get_energy(unit)
    raise TypeError(f"Unknown type: {type(ket)=}")
