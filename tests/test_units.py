# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Test the conversion from AtomicUnits used in the backend and the units input and output to the user."""

import numpy as np
import pairinteraction.real as pi
from pairinteraction.units import AtomicUnits, UnitConverterScalar, ureg


def test_magnetic() -> None:
    """Test magnetic units."""
    ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)

    mu = ket.get_matrix_element(ket, "magnetic_dipole", q=0)
    mu = mu.to("bohr_magneton")
    lande_factor = 2.002319304363
    assert np.isclose(mu.magnitude, -1 / 2 * lande_factor)

    # check magnetic field conversion is correct
    magnetic_field_au = UnitConverterScalar.user_to_au(1, "gauss", "magnetic_field")
    magnetic_field_pint = ureg.Quantity(1, "gauss").to("T", "Gaussian")
    assert np.isclose(magnetic_field_au, magnetic_field_pint.to_base_units().magnitude)

    # such that mu * magnetic_field is of dimension energy
    zeeman_energy = -mu * magnetic_field_pint
    assert zeeman_energy.dimensionality == AtomicUnits["energy"].dimensionality

    # check against constructed Hamiltonian
    basis = pi.BasisAtom("Rb", n=(1, 1), additional_kets=[ket])
    system = pi.SystemAtom(basis)
    system.set_diamagnetism_enabled(False).set_magnetic_field([0, 0, magnetic_field_pint])
    zeeman_energy_from_hamiltonian = system.get_hamiltonian("MHz")[0, 0] - ket.get_energy("MHz")
    assert np.isclose(zeeman_energy_from_hamiltonian, zeeman_energy.to("MHz", "spectroscopy").magnitude)


def test_electric_dipole() -> None:
    """Test electric dipole units."""
    ket_a = pi.KetAtom("Rb", n=60, l=0, m=0.5)
    ket_b = pi.KetAtom("Rb", n=61, l=0, m=0.5)
    ket_c = pi.KetAtom("Rb", n=60, l=1, j=3 / 2, m=0.5)

    dipole_a_c = ket_a.get_matrix_element(ket_c, "electric_dipole", q=0)
    dipole_b_c = ket_b.get_matrix_element(ket_c, "electric_dipole", q=0)

    kappa = ureg.Quantity(1 / (4 * np.pi), "1 / epsilon_0")
    c3 = kappa * dipole_a_c * dipole_b_c
    c3 = c3 * ureg.Quantity(1, "GHz") / ureg.Quantity(1, "GHz").to("hartree", "spectroscopy")
    c3 = c3.to("GHz micrometer^3")

    assert c3.magnitude == UnitConverterScalar.pint_to_user(kappa * dipole_a_c * dipole_b_c, "c3", "GHz micrometer^3")

    distance = ureg.Quantity(10, "micrometer")
    basis = pi.BasisAtom("Rb", additional_kets=[ket_a, ket_b, ket_c])
    system = pi.SystemAtom(basis)
    basis_pair = pi.BasisPair([system, system])
    system_pair = pi.SystemPair(basis_pair)
    system_pair.set_interaction_order(3)
    system_pair.set_distance(distance)
    system.get_hamiltonian()

    ket_ab_idx = np.argmax(basis_pair.get_overlaps([ket_a, ket_b]))
    ket_cc_idx = np.argmax(basis_pair.get_overlaps([ket_c, ket_c]))
    hamiltonian = system_pair.get_hamiltonian("GHz") * distance.to("micrometer").magnitude ** 3  # GHz * micrometer^3

    assert np.isclose(-2 * c3.magnitude, hamiltonian[ket_ab_idx, ket_cc_idx])
    assert np.isclose(-2 * c3.magnitude, 5.73507543166919)
