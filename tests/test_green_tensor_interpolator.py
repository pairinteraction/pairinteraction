# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from pairinteraction.units import ureg

if TYPE_CHECKING:
    from .utils import PairinteractionModule


@pytest.mark.parametrize("distance_mum", [1, 2, 11])
def test_static_green_tensor_interpolator(pi_module: PairinteractionModule, distance_mum: float) -> None:
    """Test calculating a pair potential using a user-defined static green tensor interpolator."""
    # Create a single-atom system
    basis = pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 2))
    system = pi_module.SystemAtom(basis)

    # Create two-atom basis
    ket = pi_module.KetAtom("Rb", n=60, l=0, m=0.5)
    delta_energy = 3  # GHz
    min_energy = 2 * ket.get_energy(unit="GHz") - delta_energy
    max_energy = 2 * ket.get_energy(unit="GHz") + delta_energy
    basis_pair = pi_module.BasisPair([system, system], energy=(min_energy, max_energy), energy_unit="GHz", m=(1, 1))

    # Create a system using a user-defined green tensor interpolator for dipole-dipole interaction
    gti = pi_module.GreenTensorInterpolator()
    omega = ureg.Quantity(1, "Hz")
    distance_au = ureg.Quantity(distance_mum, "micrometer").to("bohr").m

    tensor = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / (distance_au**3)
    tensor_unit = "hartree / (e * bohr)^2"
    gti.set_constant(1, 1, tensor, tensor_unit, from_coordinates="cartesian")

    tensor_spherical = gti.get(1, 1, omega, unit=tensor_unit, scaled_green_tensor=True, coordinates="spherical")
    tensor_spherical_ref = np.array([[1, 0, 0], [0, -2, 0], [0, 0, 1]]) / distance_au**3
    np.testing.assert_allclose(tensor_spherical, tensor_spherical_ref)

    tensor_spherical = gti.get(
        1, 1, ureg.Quantity(2.5, "GHz"), unit=tensor_unit, scaled_green_tensor=True, coordinates="spherical"
    )
    np.testing.assert_allclose(tensor_spherical, tensor_spherical_ref)

    system_pairs = pi_module.SystemPair(basis_pair).set_green_tensor_interpolator(gti)

    # Create a reference system using the build in dipole-dipole interaction
    system_pairs_reference = (
        pi_module.SystemPair(basis_pair).set_interaction_order(3).set_distance(distance_mum, unit="micrometer")
    )

    # Check that the two systems are equivalent
    np.testing.assert_allclose(
        system_pairs.get_hamiltonian(unit="GHz").data, system_pairs_reference.get_hamiltonian(unit="GHz").data
    )


@pytest.mark.parametrize("distance_mum", [1, 2, 11])
def test_omega_dependent_green_tensor_interpolator(pi_module: PairinteractionModule, distance_mum: float) -> None:
    """Test the interpolation for different values of omega."""
    # Define an simple linear omega-dependent green tensor interpolator
    # note that at least four entries are needed for the applied spline interpolation.
    gti = pi_module.GreenTensorInterpolator()
    omegas = np.linspace(1, 5, 20)  # GHz
    tensors = [
        np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) * i / (omega**2 * distance_mum**3)
        for (i, omega) in enumerate(omegas, start=1)
    ]
    tensor_unit = "1 / micrometer"
    gti.set_list(1, 1, tensors, omegas, tensors_unit=tensor_unit, omegas_unit="GHz", from_coordinates="cartesian")

    # Check the interpolation
    tensors_spherical = [
        np.array([[1, 0, 0], [0, -2, 0], [0, 0, 1]]) * i / (omega**2 * distance_mum**3)
        for (i, omega) in enumerate(omegas, start=1)
    ]

    for idx in range(3, len(omegas) - 3):  # the interpolation near the edges is bad, so we only check the middle
        tensor = gti.get(1, 1, omega=omegas[idx], omega_unit="GHz", unit=tensor_unit, coordinates="spherical")
        np.testing.assert_allclose(tensor, tensors_spherical[idx])

    for ind in range(3, len(omegas) - 5):
        ind1, ind2 = ind, ind + 1
        omega = (omegas[ind1] + omegas[ind2]) / 2
        reference_tensor = (tensors_spherical[ind1] + tensors_spherical[ind2]) / 2

        tensor = gti.get(1, 1, omega=omega, omega_unit="GHz", unit=tensor_unit, coordinates="spherical")
        np.testing.assert_allclose(tensor, reference_tensor, rtol=2e-2)
