# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .utils import PairinteractionModule


def test_energy_range(pi_module: PairinteractionModule) -> None:
    """Test restricting the energy range in the diagonalization."""
    ket = pi_module.KetAtom("Rb", n=60, l=0, m=0.5)
    pair_energy = 2 * ket.get_energy(unit="GHz")
    distances = np.linspace(1, 4, 10)

    # Create a single-atom system
    basis = pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 1))
    system = pi_module.SystemAtom(basis)

    # Create two-atom basis
    basis_pair = pi_module.BasisPair(
        [system, system], energy=(pair_energy - 10, pair_energy + 10), energy_unit="GHz", m=(1, 1)
    )

    # Diagonalize the systems for different distances in parallel and get all eigenenergies
    system_pairs = [pi_module.SystemPair(basis_pair).set_distance(d, unit="micrometer") for d in distances]
    pi_module.diagonalize(system_pairs, diagonalizer="eigen", sort_by_energy=True)
    eigenenergies_all = [system.get_eigenenergies(unit="GHz") for system in system_pairs]

    # Diagonalize the systems for different distances in parallel and get only the eigenenergies in an energy range
    system_pairs = [pi_module.SystemPair(basis_pair).set_distance(d, unit="micrometer") for d in distances]
    pi_module.diagonalize(
        system_pairs,
        diagonalizer="eigen",
        sort_by_energy=True,
        energy_range=(pair_energy - 5, pair_energy + 5),
        energy_range_unit="GHz",
    )
    eigenenergies_restricted = [system.get_eigenenergies(unit="GHz") for system in system_pairs]

    # Check the result
    eigenenergies_all_restricted = [
        eigenenergies[(eigenenergies < pair_energy + 5) & (eigenenergies > pair_energy - 5)]
        for eigenenergies in eigenenergies_all
    ]
    for e1, e2 in zip(eigenenergies_restricted, eigenenergies_all_restricted):
        np.testing.assert_allclose(e1, e2)
