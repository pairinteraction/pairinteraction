# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Test user-defined green tensors."""

import numpy as np
import pairinteraction.real as pi
from pairinteraction._wrapped import GreenTensorReal

from tests.constants import HARTREE_IN_GHZ, UM_IN_ATOMIC_UNITS


def test_static_green_tensor() -> None:
    """Test calculating a pair potential using a user-defined static green tensor."""
    distance = 2  # um

    # Create a single-atom system
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2))
    system = pi.SystemAtom(basis)

    # Create two-atom basis
    ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
    delta_energy = 3  # GHz
    min_energy = 2 * ket.get_energy(unit="GHz") - delta_energy
    max_energy = 2 * ket.get_energy(unit="GHz") + delta_energy
    basis_pair = pi.BasisPair([system, system], energy=(min_energy, max_energy), energy_unit="GHz", m=(1, 1))

    # Create a system using a user-defined green tensor for dipole-dipole interaction
    gt = GreenTensorReal()
    gt.set_entries(
        1,
        1,
        np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, -2],
            ]
        )
        * 1
        / (distance * UM_IN_ATOMIC_UNITS) ** 3,
    )
    system_pairs = pi.SystemPair(basis_pair).set_green_tensor(gt)

    # Create a reference system using the build in dipole-dipole interaction
    system_pairs_reference = (
        pi.SystemPair(basis_pair).set_interaction_order(3).set_distance(distance, unit="micrometer")
    )

    # Check that the two systems are equivalent
    np.testing.assert_allclose(
        system_pairs.get_hamiltonian(unit="GHz").data, system_pairs_reference.get_hamiltonian(unit="GHz").data
    )


def test_omega_dependent_green_tensor() -> None:
    """Test the interpolation for different values of omega."""
    distance = 2  # um

    # Define an omega-dependent green tensor,
    # note that at least four entries are needed for the applied spline interpolation.
    gt = GreenTensorReal()
    gt.set_entries(
        1,
        1,
        np.array(
            [
                [
                    [1, 0, 0],
                    [0, 1, 0],
                    [0, 0, -2],
                ],
                [
                    [2, 0, 0],
                    [0, 2, 0],
                    [0, 0, -4],
                ],
                [
                    [3, 0, 0],
                    [0, 3, 0],
                    [0, 0, -6],
                ],
                [
                    [4, 0, 0],
                    [0, 4, 0],
                    [0, 0, -8],
                ],
            ]
        )
        * 1
        / (distance * UM_IN_ATOMIC_UNITS) ** 3,
        [1 / HARTREE_IN_GHZ, 2 / HARTREE_IN_GHZ, 3 / HARTREE_IN_GHZ, 4 / HARTREE_IN_GHZ],
    )

    # Check the interpolation
    entries = gt.get_entries(1, 1)
    rows = [entry.row() for entry in entries]
    cols = [entry.col() for entry in entries]
    assert rows == [0, 1, 2]
    assert cols == [0, 1, 2]

    omegas = [2, 3]  # the interpolation near the edges of the range is bad, so we only check the middle
    references = [2, 3]
    for omega, reference in zip(omegas, references):
        val00 = entries[0].val(omega / HARTREE_IN_GHZ) * (distance * UM_IN_ATOMIC_UNITS) ** 3
        assert abs(val00 - reference) / reference < 0.02

    val00 = entries[0].val(2.5 / HARTREE_IN_GHZ) * (distance * UM_IN_ATOMIC_UNITS) ** 3
    assert abs(val00 - 2.5) / 2.5 < 0.02
