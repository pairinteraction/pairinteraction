# SPDX-FileCopyrightText: 2026 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .utils import PairinteractionModule


def test_state_pair_label(pi_module: PairinteractionModule) -> None:
    """Test the label of a pair eigenstate."""
    basis_atom = pi_module.BasisAtom("Rb", n=(60, 60), l=(0, 1))
    system_atom = pi_module.SystemAtom(basis_atom).set_electric_field([0, 0, 5], unit="V/cm").diagonalize()
    basis_pair = pi_module.BasisPair([system_atom, system_atom], m=(1, 1))
    system_pair = pi_module.SystemPair(basis_pair).diagonalize()

    state = system_pair.get_eigenbasis().get_state(0)
    expected_label = (
        "0.77 |Rb:60,S_1/2,1/2; Rb:60,S_1/2,1/2⟩ - 0.34 |Rb:60,P_3/2,1/2; Rb:60,S_1/2,1/2⟩ "
        "- 0.34 |Rb:60,S_1/2,1/2; Rb:60,P_3/2,1/2⟩ + ..."
    )

    assert state.get_label() == expected_label
