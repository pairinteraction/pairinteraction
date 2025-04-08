# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Test diamagnetic interactions."""

import numpy as np

import pairinteraction.complex as pi


def test_diamagnetism() -> None:
    """Test calculating diamagnetic interactions."""
    # Create a basis
    ket = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    basis = pi.BasisAtom("Rb", n=(ket.n - 2, ket.n + 2), l=(0, ket.l + 2))
    print(f"Number of basis states: {basis.number_of_states}")

    # Create system for a magnetic field of 1000 G
    bfield = 1000
    system = pi.SystemAtom(basis).set_magnetic_field([0, 0, bfield], unit="G").set_diamagnetism_enabled(True)

    # Diagonalize the system
    system = system.diagonalize(diagonalizer="eigen", sort_by_energy=True)

    # Get eigenenergies and the overlap with |ket>
    overlaps = system.basis.get_overlaps(ket)
    eigenenergies = system.get_eigenenergies(unit="GHz") - ket.get_energy(unit="GHz")

    # Get the overlap and eigenenergy corresponding to |ket>
    idx = np.argmax(overlaps)
    overlap = overlaps[idx]
    eigenenergy = eigenenergies[idx]
    print(
        f"The state |{ket}> in a field of {bfield} G has an energy of {eigenenergy:.3f} GHz "
        f"and overlap of {overlap:.2%} with the unperturbed state."
    )

    # Compare to reference data
    assert np.isclose(overlap, system.basis.get_corresponding_state(ket).get_overlap(ket))
    assert np.isclose(overlap, 0.9823116102876408, atol=1e-10)
    assert np.isclose(eigenenergy, 4.113262772909366)


def test_diamagnetism_angle_dependence() -> None:
    """Test calculating diamagnetic interactions for differently aligned magnetic fields."""
    # Create a basis
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2))

    # Create systems for fields in different directions
    bfield = 1000
    system_x = pi.SystemAtom(basis).set_magnetic_field([bfield, 0, 0], unit="G").set_diamagnetism_enabled(True)
    system_y = pi.SystemAtom(basis).set_magnetic_field([0, bfield, 0], unit="G").set_diamagnetism_enabled(True)
    system_z = pi.SystemAtom(basis).set_magnetic_field([0, 0, bfield], unit="G").set_diamagnetism_enabled(True)

    # Diagonalize the systems in parallel
    pi.diagonalize([system_x, system_y, system_z], diagonalizer="eigen", sort_by_energy=True)

    # Ensure that all eigenenergies are the same
    np.testing.assert_allclose(system_x.get_eigenenergies(unit="GHz"), system_y.get_eigenenergies(unit="GHz"))
    np.testing.assert_allclose(system_x.get_eigenenergies(unit="GHz"), system_z.get_eigenenergies(unit="GHz"))
