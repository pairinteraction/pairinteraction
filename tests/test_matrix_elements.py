"""Test the calculation of matrix elements."""

import numpy as np

import pairinteraction.backend.double as pi


def test_calculate_energy(database_dir: str, download_missing: bool) -> None:
    """Test calculating energies of ket states."""
    database = pi.Database(download_missing, True, database_dir)

    # Energy of unperturbed state
    ket = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5, database=database)
    energy_unperturbed = pi.calculate_energy(ket, unit="GHz")

    assert np.isclose(energy_unperturbed, ket.get_energy("GHz"))

    # Energy of Stark shifted state
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2), m=(0.5, 0.5), database=database)

    system = pi.SystemAtom(basis).set_electric_field([0, 0, 1], unit="V/cm").diagonalize(diagonalizer="Eigen")

    energy_perturbed = pi.calculate_energy(ket, system, unit="GHz")

    shift = energy_perturbed - energy_unperturbed
    print(f"Energy shift: {shift} GHz")

    assert shift < 0


def test_calculate_electric_dipole_matrix_element(database_dir: str, download_missing: bool) -> None:
    """Test calculating dipole matrix elements."""
    database = pi.Database(download_missing, True, database_dir)

    # The dipole element between dipole-coupled states should be non-zero
    ket_initial = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5, database=database)
    ket_final = pi.KetAtom("Rb", n=60, l=1, j=0.5, m=0.5, database=database)

    dipole = database.get_matrix_element(ket_initial, ket_final, "ELECTRIC_DIPOLE", 0, unit="e a0")
    assert dipole > 0

    # The dipole element between the same, unperturbed state should be zero
    dipole = database.get_matrix_element(ket_initial, ket_initial, "ELECTRIC_DIPOLE", 0, unit="e a0")
    assert np.isclose(dipole, 0)

    # Stark effect induces a permanent dipole moment
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2), database=database)

    system = (
        pi.SystemAtom(basis)
        .set_electric_field([1, 0, 1], unit="V/cm")
        .diagonalize(diagonalizer="Eigen", sort_by_energy=True)
    )

    dipole_zero = database.get_matrix_element(
        ket_initial, ket_initial, "ELECTRIC_DIPOLE", 0, system=system, unit="e a0"
    )
    dipole_plus = database.get_matrix_element(
        ket_initial, ket_initial, "ELECTRIC_DIPOLE", 1, system=system, unit="e a0"
    )
    dipole_minus = database.get_matrix_element(
        ket_initial, ket_initial, "ELECTRIC_DIPOLE", -1, system=system, unit="e a0"
    )

    dipole_z = dipole_zero
    dipole_x = (dipole_minus - dipole_plus) / np.sqrt(2)
    dipole_y = 1j * (dipole_plus + dipole_minus) / np.sqrt(2)

    assert dipole_z > 0
    assert np.isclose(dipole_z, dipole_x)
    assert np.isclose(dipole_y, 0)
