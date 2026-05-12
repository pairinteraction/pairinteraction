# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from pairinteraction_gui.config.base_config import BaseConfig
from pairinteraction_gui.config.basis_config import BasisConfig, BasisConfigC6, BasisConfigOneAtom, BasisConfigTwoAtoms
from pairinteraction_gui.config.calculation_config import CalculationConfig
from pairinteraction_gui.config.ket_config import (
    KetConfig,
    KetConfigC6,
    KetConfigLifetimes,
    KetConfigOneAtom,
    KetConfigTwoAtoms,
)
from pairinteraction_gui.config.system_config import (
    SystemConfig,
    SystemConfigC6,
    SystemConfigOneAtom,
    SystemConfigTwoAtoms,
)

# The naming scheme for the config classes follows the pattern:
# <ConfigType>Config<PageType>, e.g. the KetConfig for the OneAtomPage is called KetConfigOneAtom.

__all__ = [
    "BaseConfig",
    "BasisConfig",
    "BasisConfigC6",
    "BasisConfigOneAtom",
    "BasisConfigTwoAtoms",
    "CalculationConfig",
    "KetConfig",
    "KetConfigC6",
    "KetConfigLifetimes",
    "KetConfigOneAtom",
    "KetConfigTwoAtoms",
    "SystemConfig",
    "SystemConfigC6",
    "SystemConfigOneAtom",
    "SystemConfigTwoAtoms",
]
