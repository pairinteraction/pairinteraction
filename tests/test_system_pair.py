# SPDX-FileCopyrightText: 2026 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

import pytest


def test_real_system_pair_distance_vector_with_y_component_raises() -> None:
    import pairinteraction.real as pi_real

    ket = pi_real.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    basis = pi_real.BasisAtom("Rb", n=(0, 0), additional_kets=[ket])
    system = pi_real.SystemAtom(basis)
    basis_pair = pi_real.BasisPair((system, system))

    with pytest.raises(ValueError, match="y-component"):
        pi_real.SystemPair(basis_pair).set_distance_vector([0, 1, 0], unit="micrometer")
