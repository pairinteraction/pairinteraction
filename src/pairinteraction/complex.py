# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import warnings

from pairinteraction.basis import BasisAtom, BasisPair
from pairinteraction.database import Database
from pairinteraction.diagonalize import diagonalize
from pairinteraction.ket import KetAtom, KetPair
from pairinteraction.perturbative import C3, C6, EffectiveSystemPair
from pairinteraction.state import StateAtom, StatePair
from pairinteraction.system import GreenTensor, SystemAtom, SystemPair

__all__ = [
    "C3",
    "C6",
    "BasisAtom",
    "BasisPair",
    "Database",
    "EffectiveSystemPair",
    "GreenTensor",
    "KetAtom",
    "KetPair",
    "StateAtom",
    "StatePair",
    "SystemAtom",
    "SystemPair",
    "diagonalize",
]


warnings.warn(
    "All classes and functions from 'pairinteraction.complex' are now available via 'pairinteraction' directly. "
    "Replace all imports like 'import pairinteraction.complex as pi' with 'import pairinteraction as pi' "
    "to avoid this warning.",
    DeprecationWarning,
    stacklevel=2,
)
