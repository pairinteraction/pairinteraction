# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

if TYPE_CHECKING:
    from .utils import PairinteractionModule


@pytest.mark.parametrize("distance_mum", [1, 2, 11])
def test_static_green_tensor(pi_module: PairinteractionModule, distance_mum: float) -> None:
    """Test calculating a pair potential using a user-defined static green tensor."""
    # Create a single-atom system
    basis = pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 2))
    system = pi_module.SystemAtom(basis)

    # Create two-atom basis
    ket = pi_module.KetAtom("Rb", n=60, l=0, m=0.5)
    delta_energy = 3  # GHz
    min_energy = 2 * ket.get_energy(unit="GHz") - delta_energy
    max_energy = 2 * ket.get_energy(unit="GHz") + delta_energy
    basis_pair = pi_module.BasisPair([system, system], energy=(min_energy, max_energy), energy_unit="GHz", m=(1, 1))

    # Create a system using a user-defined green tensor for dipole-dipole interaction
    gt = pi_module.GreenTensor()
    tensor = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / distance_mum**3
    tensor_unit = "hartree / (e^2 micrometer^3)"
    gt.set_from_cartesian(1, 1, tensor, tensor_unit)

    tensor_spherical = np.array([[1, 0, 0], [0, -2, 0], [0, 0, 1]]) / distance_mum**3
    np.testing.assert_allclose(gt.get_spherical(1, 1, unit=tensor_unit), tensor_spherical)
    np.testing.assert_allclose(gt.get_spherical(1, 1, omega=2.5, omega_unit="GHz", unit=tensor_unit), tensor_spherical)

    system_pairs = pi_module.SystemPair(basis_pair).set_green_tensor(gt)

    # Create a reference system using the build in dipole-dipole interaction
    system_pairs_reference = (
        pi_module.SystemPair(basis_pair).set_interaction_order(3).set_distance(distance_mum, unit="micrometer")
    )

    # Check that the two systems are equivalent
    np.testing.assert_allclose(
        system_pairs.get_hamiltonian(unit="GHz").data, system_pairs_reference.get_hamiltonian(unit="GHz").data
    )


@pytest.mark.parametrize("distance_mum", [1, 2, 11])
def test_omega_dependent_green_tensor(pi_module: PairinteractionModule, distance_mum: float) -> None:
    """Test the interpolation for different values of omega."""
    # Define an simple linear omega-dependent green tensor
    # note that at least four entries are needed for the applied spline interpolation.
    gt = pi_module.GreenTensor()
    omegas = [1, 2, 3, 4]  # GHz
    tensors = [np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) * i / distance_mum**3 for i in range(1, 5)]
    tensor_unit = "hartree / (e^2 micrometer^3)"
    gt.set_from_cartesian(1, 1, tensors, tensor_unit, omegas, omegas_unit="GHz")

    # Check the interpolation
    tensors_spherical = [np.array([[1, 0, 0], [0, -2, 0], [0, 0, 1]]) * i / distance_mum**3 for i in range(1, 5)]

    for idx in range(1, len(omegas) - 1):  # the interpolation near the edges is bad, so we only check the middle
        tensor = gt.get_spherical(1, 1, omega=omegas[idx], omega_unit="GHz", unit=tensor_unit)
        np.testing.assert_allclose(tensor, tensors_spherical[idx])

    tensor = gt.get_spherical(1, 1, omega=2.5, omega_unit="GHz", unit=tensor_unit)
    reference_tensor = (tensors_spherical[1] + tensors_spherical[2]) / 2
    np.testing.assert_allclose(tensor, reference_tensor, rtol=2e-2)
