# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Test atom-ion interaction."""

import numpy as np

import pairinteraction.complex as pi


def test_ion() -> None:
    """Test the calculation of energy shifts in the field on an ion."""
    # Create a basis
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2))

    # Create systems for different distances to the ion
    distances = np.linspace(1, 10, 5)
    systems_x = [
        pi.SystemAtom(basis).set_ion_charge(1, unit="e").set_distance_vector_to_ion([d, 0, 0], unit="um")
        for d in distances
    ]
    systems_y = [
        pi.SystemAtom(basis).set_ion_charge(1, unit="e").set_distance_vector_to_ion([0, d, 0], unit="um")
        for d in distances
    ]
    systems_z = [
        pi.SystemAtom(basis).set_ion_charge(1, unit="e").set_distance_vector_to_ion([0, 0, d], unit="um")
        for d in distances
    ]

    # Diagonalize the systems in parallel
    pi.diagonalize(systems_x, diagonalizer="eigen", sort_by_energy=True)
    pi.diagonalize(systems_y, diagonalizer="eigen", sort_by_energy=True)
    pi.diagonalize(systems_z, diagonalizer="eigen", sort_by_energy=True)

    # Ensure that all eigenenergies are the same
    for system_x, system_y, system_z in zip(systems_x, systems_y, systems_z):
        np.testing.assert_allclose(system_x.get_eigenenergies(unit="GHz"), system_y.get_eigenenergies(unit="GHz"))
        np.testing.assert_allclose(system_x.get_eigenenergies(unit="GHz"), system_z.get_eigenenergies(unit="GHz"))
