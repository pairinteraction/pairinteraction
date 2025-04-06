# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Test atom-ion interaction."""

import numpy as np
import pytest

import pairinteraction.complex as pi


def test_ion_z() -> None:
    """Test the calculation of energy shifts in the field on an ion positioned along z."""
    # Create a basis
    ket = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    basis = pi.BasisAtom("Rb", n=(ket.n - 2, ket.n + 2), l=(0, ket.l + 2), m=(ket.m, ket.m))

    # Create systems for different distances to the ion
    distance = 3
    system_z = (
        pi.SystemAtom(basis)
        .set_ion_interaction_order(3)
        .set_ion_charge(1, unit="e")
        .set_ion_distance_vector([0, 0, distance], unit="um")
    )

    # Diagonalize the system
    system_z = system_z.diagonalize(diagonalizer="eigen", sort_by_energy=True)

    # Ensure that values are correct for the system where the ion is closest to the atom
    eigenenergies = system_z.get_eigenenergies(unit="GHz")
    overlaps = system_z.basis.get_overlaps(ket)
    idx = np.argmax(overlaps)
    assert pytest.approx(overlaps[idx], rel=1e-6) == 0.8754278674619744  # NOSONAR
    assert pytest.approx(eigenenergies[idx] - ket.get_energy(unit="GHz"), rel=1e-12) == -0.2531503909267485  # NOSONAR


def test_ion_x() -> None:
    """Test the calculation of energy shifts in the field on an ion positioned along x."""
    # Create a basis
    ket = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    basis = pi.BasisAtom("Rb", n=(ket.n - 2, ket.n + 2), l=(0, ket.l + 2))

    # Create systems for different distances to the ion
    distance = 3
    system_x = (
        pi.SystemAtom(basis)
        .set_ion_interaction_order(3)
        .set_ion_charge(1, unit="e")
        .set_ion_distance_vector([distance, 0, 0], unit="um")
    )

    # Diagonalize the system
    system_x = system_x.diagonalize(diagonalizer="eigen", sort_by_energy=True)

    # Ensure that values are correct for the system where the ion is closest to the atom
    eigenenergies = system_x.get_eigenenergies(unit="GHz")
    overlaps = system_x.basis.get_overlaps(ket)
    idx = np.argmax(overlaps)
    # Note that we must use a large relative tolerance for the overlaps because the calculated value
    # is very sensitive on the actual method that is used by eigen to diagonalize the system (whether eigen
    # is using its own implementation, mkl on a Intel CPU, mkl on a AMD CPU, or lapack). This is because of
    # eigenstates belonging to different degenerate Zeeman sublevels.
    assert pytest.approx(overlaps[idx], rel=0.2) == 0.8754157131466024  # NOSONAR
    assert pytest.approx(eigenenergies[idx] - ket.get_energy(unit="GHz"), rel=1e-12) == -0.2531503909267485  # NOSONAR


def test_ion_angle_dependence() -> None:
    """Test the calculation of energy shifts in the field on an ion for different angles."""
    # Create a basis
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2))

    # Create systems for different distances to the ion
    distances = np.linspace(3, 10, 5)
    systems_x = [
        pi.SystemAtom(basis)
        .set_ion_interaction_order(3)
        .set_ion_charge(1, unit="e")
        .set_ion_distance_vector([d, 0, 0], unit="um")
        for d in distances
    ]
    systems_y = [
        pi.SystemAtom(basis)
        .set_ion_interaction_order(3)
        .set_ion_charge(1, unit="e")
        .set_ion_distance_vector([0, d, 0], unit="um")
        for d in distances
    ]
    systems_z = [
        pi.SystemAtom(basis)
        .set_ion_interaction_order(3)
        .set_ion_charge(1, unit="e")
        .set_ion_distance_vector([0, 0, d], unit="um")
        for d in distances
    ]

    # Diagonalize the systems in parallel
    pi.diagonalize(systems_x + systems_y + systems_z, diagonalizer="eigen", sort_by_energy=True)

    # Ensure that all eigenenergies are the same
    for system_x, system_y, system_z in zip(systems_x, systems_y, systems_z):
        np.testing.assert_allclose(system_x.get_eigenenergies(unit="GHz"), system_y.get_eigenenergies(unit="GHz"))
        np.testing.assert_allclose(system_x.get_eigenenergies(unit="GHz"), system_z.get_eigenenergies(unit="GHz"))
