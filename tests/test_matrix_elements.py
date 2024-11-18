"""
Test the calculation of matrix elements.
"""

import numpy as np

import pairinteraction.backend.double as pi
from pairinteraction import ureg


def test_calculate_energy(database_dir: str, download_missing: bool) -> None:
    """Test calculating energies of ket states."""
    database = pi.Database(download_missing, True, database_dir)
    diagonalizer = pi.DiagonalizerEigen()

    # Energy of unperturbed state
    ket = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5, database=database)
    energy_unperturbed = pi.calculate_energy(ket)

    assert np.isclose(energy_unperturbed, ket.get_energy())

    # Energy of Stark shifted state
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2), m=(0.5, 0.5), database=database)

    system = pi.SystemAtom(basis).set_electric_field([0, 0, 1], unit="V/cm").diagonalize(diagonalizer)

    energy_perturbed = pi.calculate_energy(ket, system)

    shift = (
        (energy_perturbed * ureg.hartree - energy_unperturbed * ureg.hartree)
        .to(ureg.gigahertz, "spectroscopy")
        .magnitude
    )
    print(f"Energy shift: {shift} GHz")

    assert shift < 0


def test_calculate_electric_dipole_matrix_element(database_dir: str, download_missing: bool) -> None:
    """Test calculating dipole matrix elements."""
    database = pi.Database(download_missing, True, database_dir)
    diagonalizer = pi.DiagonalizerEigen()

    # The dipole element between dipole-coupled states should be non-zero
    ket_initial = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5, database=database)
    ket_final = pi.KetAtom("Rb", n=60, l=1, j=0.5, m=0.5, database=database)

    dipole = pi.calculate_electric_dipole_matrix_element(ket_initial, ket_final, 0)
    assert dipole > 0

    # The dipole element between the same, unperturbed state should be zero
    ket_final = ket_initial

    dipole = pi.calculate_electric_dipole_matrix_element(ket_initial, ket_final, 0)
    assert np.isclose(dipole, 0)

    # Stark effect induces a permanent dipole moment
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2), database=database)

    system = pi.SystemAtom(basis).set_electric_field([1, 0, 1], unit="V/cm").diagonalize(diagonalizer)

    dipole_zero = pi.calculate_electric_dipole_matrix_element(ket_initial, ket_final, system, 0)
    dipole_plus = pi.calculate_electric_dipole_matrix_element(ket_initial, ket_final, system, 1)
    dipole_minus = pi.calculate_electric_dipole_matrix_element(ket_initial, ket_final, system, -1)

    dipole_z = dipole_zero
    dipole_x = (dipole_minus - dipole_plus) / np.sqrt(2)
    dipole_y = 1j * (dipole_plus + dipole_minus) / np.sqrt(2)

    assert dipole_z > 0
    assert np.isclose(dipole_z, dipole_x)
    assert np.isclose(dipole_y, 0)
