# SPDX-FileCopyrightText: 2026 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

import numpy as np
import pytest


def test_real_system_pair_distance_vector_with_y_component_raises() -> None:
    import pairinteraction.real as pi_real

    ket = pi_real.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    basis = pi_real.BasisAtom("Rb", n=(0, 0), additional_kets=[ket])
    system = pi_real.SystemAtom(basis)
    basis_pair = pi_real.BasisPair((system, system))

    with pytest.raises(ValueError, match="y-component"):
        pi_real.SystemPair(basis_pair).set_distance_vector([0, 1, 0], unit="micrometer")


def test_effective_system_pair_set_distance_pint() -> None:
    import pairinteraction.real as pi_real

    ket = pi_real.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    eff_system = pi_real.EffectiveSystemPair([(ket, ket)])
    distance = pi_real.ureg.Quantity(5, "micrometer")
    expected = [2.5, 0, 5 * np.cos(np.pi / 6)]

    eff_system.set_distance(distance, angle_degree=30)
    np.testing.assert_allclose(eff_system.distance_vector.to("micrometer").magnitude, expected)

    eff_system.set_distance_vector([distance * 0.5, 0, distance * np.cos(np.pi / 6)])
    np.testing.assert_allclose(eff_system.distance_vector.to("micrometer").magnitude, expected)

    eff_system.set_distance_vector(pi_real.ureg.Quantity(np.array(expected), "micrometer"))
    np.testing.assert_allclose(eff_system.distance_vector.to("micrometer").magnitude, expected)

    eff_system.set_distance_vector(expected, "micrometer")
    np.testing.assert_allclose(eff_system.distance_vector.to("micrometer").magnitude, expected)
