# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from collections.abc import Sequence
from functools import cached_property
from typing import Any, ClassVar, Literal, Union

from typing_extensions import TypeGuard

from pairinteraction import _backend
from pairinteraction.basis.basis_atom import BasisAtomComplex, BasisAtomReal
from pairinteraction.ket.ket import KetBase
from pairinteraction.ket.ket_atom import KetAtom
from pairinteraction.state.state_atom import StateAtomComplex, StateAtomReal

UnionCPPKetPair = Union[_backend.KetPairReal, _backend.KetPairComplex]
UnionTypeCPPKetPairCreator = Any

KetAtomTuple = Union[tuple["KetAtom", "KetAtom"], Sequence["KetAtom"]]
KetPairLike = Union["KetPairReal", "KetPairComplex", KetAtomTuple]


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

    _cpp: UnionCPPKetPair
    _cpp_creator: ClassVar[UnionTypeCPPKetPairCreator] = None
    _TypeBasisAtom: Union[type[BasisAtomReal], type[BasisAtomComplex]]  # should be ClassVar, but cannot be nested yet
    _TypeStateAtom: Union[type[StateAtomReal], type[StateAtomComplex]]  # should be ClassVar, but cannot be nested yet

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

    @cached_property
    def state_atoms(self) -> Union[tuple[StateAtomReal, StateAtomReal], tuple[StateAtomComplex, StateAtomComplex]]:
        """Return the state atoms of the ket pair."""
        state_atoms = []
        for atomic_state in self._cpp.get_atomic_states():
            basis = self._TypeBasisAtom._from_cpp_object(atomic_state)
            state = self._TypeStateAtom._from_basis_object(basis)
            state_atoms.append(state)
        return tuple(state_atoms)  # type: ignore [return-value]


class KetPairReal(KetPair):
    _cpp: _backend.KetPairReal
    _TypeBasisAtom = BasisAtomReal
    _TypeStateAtom = StateAtomReal
    state_atoms: tuple[StateAtomReal, StateAtomReal]


class KetPairComplex(KetPair):
    _cpp: _backend.KetPairComplex
    _TypeBasisAtom = BasisAtomComplex
    _TypeStateAtom = StateAtomComplex
    state_atoms: tuple[StateAtomComplex, StateAtomComplex]
