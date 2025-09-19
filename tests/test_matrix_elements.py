# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .utils import PairinteractionModule


def test_energy(pi_module: PairinteractionModule) -> None:
    """Test calculating energies of ket states."""
    # Energy of unperturbed state
    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    energy_unperturbed = ket.get_energy(unit="GHz")

    assert np.isclose(energy_unperturbed, ket.get_energy("GHz"))

    # Energy of Stark shifted state
    basis = pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 2), m=(0.5, 0.5))

    system = pi_module.SystemAtom(basis).set_electric_field([0, 0, 1], unit="V/cm").diagonalize(diagonalizer="eigen")

    energy_perturbed = system.get_corresponding_energy(ket, unit="GHz")

    shift = energy_perturbed - energy_unperturbed
    print(f"Energy shift: {shift} GHz")

    assert shift < 0


def test_electric_dipole_matrix_element(pi_module: PairinteractionModule) -> None:
    """Test calculating dipole matrix elements."""
    # The dipole element between dipole-coupled states should be non-zero
    ket_initial = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket_final = pi_module.KetAtom("Rb", n=60, l=1, j=0.5, m=0.5)

    dipole = ket_initial.get_matrix_element(ket_final, "electric_dipole", 0, unit="e a0")
    assert dipole.real > 0
    assert dipole.imag == 0

    # The dipole element between the same, unperturbed state should be zero
    dipole = ket_initial.get_matrix_element(ket_initial, "electric_dipole", 0, unit="e a0")
    assert np.isclose(dipole, 0)

    # Stark effect induces a permanent dipole moment
    basis = pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 2))

    system = (
        pi_module.SystemAtom(basis)
        .set_electric_field([1, 0, 1], unit="V/cm")
        .diagonalize(diagonalizer="eigen", sort_by_energy=True)
    )

    state = system.basis.get_corresponding_state(ket_initial)

    dipole_zero = state.get_matrix_element(state, "electric_dipole", 0, unit="e a0")
    dipole_plus = state.get_matrix_element(state, "electric_dipole", 1, unit="e a0")
    dipole_minus = state.get_matrix_element(state, "electric_dipole", -1, unit="e a0")

    dipole_z = dipole_zero
    dipole_x = (dipole_minus - dipole_plus) / np.sqrt(2)
    dipole_y = 1j * (dipole_plus + dipole_minus) / np.sqrt(2)

    assert dipole_z.real > 0
    assert dipole_z.imag == 0
    assert np.isclose(dipole_z, dipole_x)
    assert np.isclose(dipole_y, 0)
