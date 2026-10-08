# SPDX-FileCopyrightText: 2026 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

import logging

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


def test_get_corresponding_energy_symmetrized(caplog: pytest.LogCaptureFixture) -> None:
    import pairinteraction as pi

    ket_a = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket_b = pi.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    basis_atom = pi.BasisAtom("Rb", n=(58, 63), l=(0, 2), m=(-1.5, 1.5))
    system_atom = pi.SystemAtom(basis_atom).diagonalize(diagonalizer="eigen")
    pair_energy = system_atom.get_corresponding_energy(ket_a, "GHz") + system_atom.get_corresponding_energy(
        ket_b, "GHz"
    )
    basis_pair = pi.BasisPair(
        [system_atom, system_atom], m=(1, 1), energy=(pair_energy - 3, pair_energy + 3), energy_unit="GHz"
    )
    assert np.allclose(np.sort(basis_pair.get_overlaps((ket_a, ket_b)))[-2:], 0.5)

    # Without interaction, the ket is distributed over two degenerate states, so the energy is still unique
    system_pair = pi.SystemPair(basis_pair).diagonalize(diagonalizer="eigen")
    with caplog.at_level(logging.WARNING):
        energy = system_pair.get_corresponding_energy((ket_a, ket_b), "GHz")
    assert "Cannot find the uniquely corresponding" not in caplog.text
    assert np.isclose(energy, pair_energy, rtol=1e-14)

    # With interaction, the symmetric and antisymmetric states have different energies
    caplog.clear()
    system_pair = pi.SystemPair(basis_pair).set_distance(5, unit="micrometer").diagonalize(diagonalizer="eigen")
    with caplog.at_level(logging.WARNING):
        system_pair.get_corresponding_energy((ket_a, ket_b), "GHz")
    assert "Cannot find the uniquely corresponding energy" in caplog.text
