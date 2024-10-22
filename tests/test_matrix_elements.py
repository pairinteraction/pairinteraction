"""
Test the calculation of matrix elements.
"""

import numpy as np
from pint import UnitRegistry

from pairinteraction.backend import (
    BasisAtomCreatorDouble,
    Database,
    DiagonalizerEigenDouble,
    KetAtomCreatorDouble,
    SystemAtomDouble,
    calculate_electric_dipole_matrix_element,
    calculate_energy,
)


def test_calculate_energy(ureg: UnitRegistry, database_dir: str, download_missing: bool) -> None:
    """Test calculating energies of ket states."""
    database = Database(download_missing, True, database_dir)
    diagonalizer = DiagonalizerEigenDouble()

    # Energy of unperturbed state
    ket = (
        KetAtomCreatorDouble()
        .set_species("Rb")
        .set_quantum_number_n(60)
        .set_quantum_number_l(0)
        .set_quantum_number_j(0.5)
        .set_quantum_number_m(0.5)
        .create(database)
    )

    energy_unperturbed = calculate_energy(ket)

    assert np.isclose(energy_unperturbed, ket.get_energy())

    # Energy of Stark shifted state
    basis = (
        BasisAtomCreatorDouble()
        .set_species("Rb")
        .restrict_quantum_number_n(58, 62)
        .restrict_quantum_number_l(0, 2)
        .restrict_quantum_number_m(0.5, 0.5)
        .create(database)
    )

    field = 1 * ureg.volt / ureg.centimeter
    system = (
        SystemAtomDouble(basis).set_electric_field([0, 0, field.to_base_units().magnitude]).diagonalize(diagonalizer)
    )

    energy_perturbed = calculate_energy(ket, system)

    shift = (
        (energy_perturbed * ureg.hartree - energy_unperturbed * ureg.hartree)
        .to(ureg.gigahertz, "spectroscopy")
        .magnitude
    )
    print(f"Energy shift: {shift} GHz")

    assert shift < 0


def test_calculate_electric_dipole_matrix_element(
    ureg: UnitRegistry, database_dir: str, download_missing: bool
) -> None:
    """Test calculating dipole matrix elements."""
    database = Database(download_missing, True, database_dir)
    diagonalizer = DiagonalizerEigenDouble()

    # The dipole element between dipole-coupled states should be non-zero
    ket_initial = (
        KetAtomCreatorDouble()
        .set_species("Rb")
        .set_quantum_number_n(60)
        .set_quantum_number_l(0)
        .set_quantum_number_j(0.5)
        .set_quantum_number_m(0.5)
        .create(database)
    )
    ket_final = (
        KetAtomCreatorDouble()
        .set_species("Rb")
        .set_quantum_number_n(60)
        .set_quantum_number_l(1)
        .set_quantum_number_j(0.5)
        .set_quantum_number_m(0.5)
        .create(database)
    )

    dipole = calculate_electric_dipole_matrix_element(ket_initial, ket_final, 0)
    assert dipole > 0

    # The dipole element between the same, unperturbed state should be zero
    ket_final = ket_initial

    dipole = calculate_electric_dipole_matrix_element(ket_initial, ket_final, 0)
    assert np.isclose(dipole, 0)

    # Stark effect induces a permanent dipole moment
    basis = (
        BasisAtomCreatorDouble()
        .set_species("Rb")
        .restrict_quantum_number_n(58, 62)
        .restrict_quantum_number_l(0, 2)
        .create(database)
    )

    field = 1 * ureg.volt / ureg.centimeter
    system = (
        SystemAtomDouble(basis)
        .set_electric_field([field.to_base_units().magnitude, 0, field.to_base_units().magnitude])
        .diagonalize(diagonalizer)
    )

    dipole_zero = calculate_electric_dipole_matrix_element(ket_initial, ket_final, system, 0)
    dipole_plus = calculate_electric_dipole_matrix_element(ket_initial, ket_final, system, 1)
    dipole_minus = calculate_electric_dipole_matrix_element(ket_initial, ket_final, system, -1)

    dipole_z = dipole_zero
    dipole_x = (dipole_minus - dipole_plus) / np.sqrt(2)
    dipole_y = 1j * (dipole_plus + dipole_minus) / np.sqrt(2)

    assert dipole_z > 0
    assert np.isclose(dipole_z, dipole_x)
    assert np.isclose(dipole_y, 0)
