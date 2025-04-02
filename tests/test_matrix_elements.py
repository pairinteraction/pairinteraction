# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Test the calculation of matrix elements."""

import numpy as np

import pairinteraction.real as pi


def test_energy() -> None:
    """Test calculating energies of ket states."""
    # Energy of unperturbed state
    ket = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    energy_unperturbed = ket.get_energy(unit="GHz")

    assert np.isclose(energy_unperturbed, ket.get_energy("GHz"))

    # Energy of Stark shifted state
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2), m=(0.5, 0.5))

    system = pi.SystemAtom(basis).set_electric_field([0, 0, 1], unit="V/cm").diagonalize(diagonalizer="eigen")

    energy_perturbed = system.get_corresponding_energy(ket, unit="GHz")

    shift = energy_perturbed - energy_unperturbed
    print(f"Energy shift: {shift} GHz")

    assert shift < 0


def test_electric_dipole_matrix_element() -> None:
    """Test calculating dipole matrix elements."""
    # The dipole element between dipole-coupled states should be non-zero
    ket_initial = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket_final = pi.KetAtom("Rb", n=60, l=1, j=0.5, m=0.5)

    dipole = ket_initial.get_matrix_element(ket_final, "electric_dipole", 0, unit="e a0")
    assert dipole > 0

    # The dipole element between the same, unperturbed state should be zero
    dipole = ket_initial.get_matrix_element(ket_initial, "electric_dipole", 0, unit="e a0")
    assert np.isclose(dipole, 0)

    # Stark effect induces a permanent dipole moment
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2))

    system = (
        pi.SystemAtom(basis)
        .set_electric_field([1, 0, 1], unit="V/cm")
        .diagonalize(diagonalizer="eigen", sort_by_energy=True)
    )

    state = system.basis.get_corresponding_state(ket_initial)

    dipole_zero = state.get_matrix_elements(state, "electric_dipole", 0, unit="e a0").todense()[0, 0]
    dipole_plus = state.get_matrix_elements(state, "electric_dipole", 1, unit="e a0").todense()[0, 0]
    dipole_minus = state.get_matrix_elements(state, "electric_dipole", -1, unit="e a0").todense()[0, 0]

    dipole_z = dipole_zero
    dipole_x = (dipole_minus - dipole_plus) / np.sqrt(2)
    dipole_y = 1j * (dipole_plus + dipole_minus) / np.sqrt(2)

    assert dipole_z > 0
    assert np.isclose(dipole_z, dipole_x)
    assert np.isclose(dipole_y, 0)
