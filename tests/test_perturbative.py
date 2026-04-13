# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import pytest
from pairinteraction import ureg
from pairinteraction.perturbative.perturbation_theory import calculate_perturbative_hamiltonian
from scipy import sparse

from .utils import no_log_propagation

if TYPE_CHECKING:
    from pairinteraction import SystemPair
    from scipy.sparse import csr_matrix

    from .utils import PairinteractionModule


def _check_sparse_matrices_equal(matrix_a: csr_matrix, matrix_b: csr_matrix) -> bool:
    """Check for equality of sparse matrices efficiently.

    This functions compares two sparse matrices in compressed sparse row format on their equality.

    Args:
        matrix_a: A sparse matrix in csr format.
        matrix_b: A sparse matrix in csr format.

    Returns:
        bool: True if matrices are equal, False if not.

    """
    matrix_a.sort_indices()
    matrix_b.sort_indices()
    if not (
        matrix_a.format == "csr"
        and matrix_b.format == "csr"
        and len(matrix_a.indices) == len(matrix_b.indices)
        and len(matrix_a.indptr) == len(matrix_b.indptr)
        and len(matrix_a.data) == len(matrix_b.data)
    ):
        return False

    return bool(
        np.all(matrix_a.indices == matrix_b.indices)
        and np.all(matrix_a.indptr == matrix_b.indptr)
        and np.allclose(matrix_a.data, matrix_b.data, rtol=0, atol=1e-14)
    )


@pytest.fixture
def system_pair_sample(pi_module: PairinteractionModule) -> SystemPair:
    basis = pi_module.BasisAtom(
        species="Rb",
        n=(59, 63),
        l=(0, 1),
        m=(-1.5, 1.5),
    )
    system = pi_module.SystemAtom(basis=basis)
    system.set_diamagnetism_enabled(False)
    system.set_magnetic_field([0, 0, 1e-3], "gauss")
    if not system.is_diagonal:
        pi_module.diagonalize([system], diagonalizer="eigen", sort_by_energy=False)
    basis_pair = pi_module.BasisPair([system, system])
    system_pair = pi_module.SystemPair(basis_pair)
    theta = 0
    r = 12
    system_pair.set_distance_vector(r * np.array([np.sin(theta), 0, np.cos(theta)]), "micrometer")
    system_pair.set_interaction_order(3)
    return system_pair


def test_perturbative_calculation1(caplog: pytest.LogCaptureFixture) -> None:
    """Test of mathematic functionality."""
    hamiltonian = sparse.csr_matrix(
        [[0, 1, 1, 0, 2], [1, 1, 0, 1, 3], [1, 0, 10, 1, 0], [0, 1, 1, 11, 1], [2, 3, 0, 1, 12]]
    )
    model_space_indices = [0, 1]

    # Order 0
    h_eff_dict, eig_perturb = calculate_perturbative_hamiltonian(hamiltonian, model_space_indices, perturbation_order=0)
    h_eff = sum(h for h in h_eff_dict.values())
    assert np.any(h_eff == np.array([[0, 0], [0, 1]]))
    assert _check_sparse_matrices_equal(eig_perturb, sparse.csr_matrix(sparse.eye(2, 5, k=0)))

    # Order 1
    h_eff_dict, eig_perturb = calculate_perturbative_hamiltonian(hamiltonian, model_space_indices, perturbation_order=1)
    h_eff = sum(h for h in h_eff_dict.values())
    assert np.any(h_eff == np.array([[0, 1], [1, 1]]))
    assert _check_sparse_matrices_equal(eig_perturb, sparse.csr_matrix([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]]))

    # Order 2
    h_eff_dict, eig_perturb = calculate_perturbative_hamiltonian(hamiltonian, model_space_indices, perturbation_order=2)
    h_eff = sum(h for h in h_eff_dict.values())
    a_00 = 0 + 1 * 1 / (0 - 10) + 2 * 2 / (0 - 12)
    a_01 = 1 + 2 * 3 / (0 - 12)
    a_10 = 1 + 3 * 2 / (1 - 12)
    a_11 = 1 + 1 * 1 / (1 - 11) + 3 * 3 / (1 - 12)
    hamiltonian_new = np.array([[a_00, (a_01 + a_10) / 2], [(a_01 + a_10) / 2, a_11]])
    assert np.any(h_eff == hamiltonian_new)

    v0 = sparse.eye(1, 5, k=0) + 1 / (0 - 10) * sparse.eye(1, 5, k=2) + 2 / (0 - 12) * sparse.eye(1, 5, k=4)
    v1 = sparse.eye(1, 5, k=1) + 1 / (1 - 11) * sparse.eye(1, 5, k=3) + 3 / (1 - 12) * sparse.eye(1, 5, k=4)
    assert _check_sparse_matrices_equal(eig_perturb, sparse.csr_matrix(sparse.vstack([v0, v1])))

    # Order 3
    with caplog.at_level(logging.ERROR):
        h_eff_dict, eig_perturb = calculate_perturbative_hamiltonian(
            hamiltonian, model_space_indices, perturbation_order=3
        )
    h_eff = sum(h for h in h_eff_dict.values())
    a_00 -= 2 * 3 * 1 / ((0 - 12) * (1 - 12))
    a_01 += (
        1 * 1 * 1 / ((1 - 10) * (1 - 11))
        + 2 * 1 * 1 / ((1 - 11) * (1 - 12))
        - 1 * 1 * 1 / ((0 - 10) * (1 - 10))
        - 2 * 2 * 1 / ((1 - 12) * (0 - 12))
    )
    a_10 += (
        1 * 1 * 1 / ((0 - 10) * (0 - 11))
        + 2 * 1 * 1 / ((0 - 11) * (0 - 12))
        - 1 * 1 * 1 / ((0 - 11) * (1 - 11))
        - 3 * 3 * 1 / ((0 - 12) * (1 - 12))
    )
    a_11 += 1 * 1 * 3 / ((1 - 11) * (1 - 12)) + 3 * 1 * 1 / ((1 - 11) * (1 - 12)) - 3 * 2 * 1 / ((1 - 12) * (0 - 12))
    hamiltonian_new = np.array([[a_00, (a_01 + a_10) / 2], [(a_01 + a_10) / 2, a_11]])
    assert np.any(h_eff == hamiltonian_new)

    v0 = v0 + (
        -0.5 * (1 * 1 / (0 - 10) ** 2 + 2 * 2 / (0 - 12) ** 2) * sparse.eye(1, 5, k=0)
        + (1 * 1 / ((0 - 10) * (0 - 11)) + 2 * 1 / ((0 - 11) * (0 - 12))) * sparse.eye(1, 5, k=3)
        + 1 * 1 / ((0 - 11) * (0 - 1)) * sparse.eye(1, 5, k=3)
        + 3 * 1 / ((0 - 1) * (0 - 12)) * sparse.eye(1, 5, k=4)
        + 3 * 2 / ((0 - 1) * (0 - 12)) * sparse.eye(1, 5, k=1)
    )
    v1 = v1 + (
        -0.5 * (1 * 1 / (1 - 11) ** 2 + 3 * 3 / (1 - 12) ** 2) * sparse.eye(1, 5, k=1)
        + 1 * 1 / ((1 - 10) * (1 - 11)) * sparse.eye(1, 5, k=2)
        + 1 * 3 / ((1 - 11) * (1 - 12)) * sparse.eye(1, 5, k=3)
        + 1 * 1 / ((1 - 11) * (1 - 12)) * sparse.eye(1, 5, k=4)
        + 1 * 1 / ((1 - 0) * (1 - 10)) * sparse.eye(1, 5, k=2)
        + 2 * 1 / ((1 - 0) * (1 - 12)) * sparse.eye(1, 5, k=4)
        + 3 * 2 / ((1 - 0) * (1 - 12)) * sparse.eye(1, 5, k=0)
    )
    assert _check_sparse_matrices_equal(eig_perturb, sparse.csr_matrix(sparse.vstack([v0, v1])))


def test_c3_with_sample_system(pi_module: PairinteractionModule, system_pair_sample: SystemPair) -> None:
    """Test whether the C3 coefficient with a given system is calculated correctly."""
    ket1 = pi_module.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    ket2 = pi_module.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5)
    c3_obj = pi_module.C3(ket1, ket2)
    c3_obj._distance_vector = None  # avoid warning due when setting system pair
    c3_obj.set_system_pair(system_pair_sample)

    c3 = c3_obj.get(unit="planck_constant * gigahertz * micrometer^3")
    assert np.isclose(-0.5 * c3, 3.1515)


def test_c3(pi_module: PairinteractionModule) -> None:
    """Test whether the C3 coefficient with automatically constructed system is calculated correctly."""
    ket1 = pi_module.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    ket2 = pi_module.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5)
    c3_obj = pi_module.C3(ket1, ket2)

    c3_obj.set_electric_field([0, 0, 0], "volt/cm")
    c3_obj.set_magnetic_field([0, 0, 10], "gauss")

    c3 = c3_obj.get(unit="planck_constant * gigahertz * micrometer^3")
    assert np.isclose(-0.5 * c3, 3.2188)


def test_c6_with_sample_system(pi_module: PairinteractionModule, system_pair_sample: SystemPair) -> None:
    """Test whether the C6 coefficient with a given system is calculated correctly."""
    ket = pi_module.KetAtom(species="Rb", n=61, l=0, j=0.5, m=0.5)
    c6_obj = pi_module.C6(ket, ket)
    c6_obj._distance_vector = None  # avoid warning due when setting system pair
    c6_obj.set_system_pair(system_pair_sample)

    c6 = c6_obj.get(unit="planck_constant * gigahertz * micrometer^6")
    assert np.isclose(c6, -167.880)


def test_c6(pi_module: PairinteractionModule) -> None:
    """Test whether the C6 coefficient with automatically constructed system is calculated correctly."""
    ket = pi_module.KetAtom(species="Rb", n=61, l=0, j=0.5, m=0.5)
    c6_obj = pi_module.C6(ket, ket)

    c6_obj.set_electric_field([0, 0, 0], "volt/cm")
    c6_obj.set_magnetic_field([0, 0, 10], "gauss")

    c6 = c6_obj.get(unit="planck_constant * gigahertz * micrometer^6")
    assert np.isclose(c6, -169.149)


def test_exact_resonance_detection(
    pi_module: PairinteractionModule, system_pair_sample: SystemPair, capsys: pytest.CaptureFixture[str]
) -> None:
    """Test whether resonance with infinite admixture is correctly detected."""
    ket1 = pi_module.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    ket2 = pi_module.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5)
    eff_system = pi_module.EffectiveSystemPair([(ket1, ket2)])
    eff_system.set_system_pair(system_pair_sample)

    # workaround to test for errors, without showing them in the std output
    with no_log_propagation("pairinteraction"), np.errstate(invalid="ignore"):
        eff_system.get_effective_hamiltonian()
    captured = capsys.readouterr()
    assert "gets a large dressing (inf overlap)" in captured.err
    assert "|~Rb:61,P_3/2,1/2; Rb:61,S_1/2,1/2⟩ has admixture inf" in captured.err


def test_near_resonance_detection(pi_module: PairinteractionModule, capsys: pytest.CaptureFixture[str]) -> None:
    """Test whether a near resonance is correctly detected."""
    ket1 = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket2 = pi_module.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    eff_system = pi_module.EffectiveSystemPair([(ket1, ket2), (ket2, ket1)])
    eff_system.set_magnetic_field([0, 0, 245], "gauss")
    eff_system.set_distance(10, 35.1, "micrometer")
    eff_system.create_basis_pair(number_of_kets=100)

    # workaround to test for errors, without showing them in the std output
    with no_log_propagation("pairinteraction"):
        eff_system.get_effective_hamiltonian()
        eff_system.check_for_resonances(0.01)
    captured = capsys.readouterr()
    assert "The most perturbing states are" in captured.err
    assert "Rb:60,P_3/2,1/2; Rb:60,P_3/2,3/2" in captured.err


def test_set_basis_atoms(pi_module: PairinteractionModule) -> None:
    ket = pi_module.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    basis0 = pi_module.BasisAtom("Rb", n=(58, 64), l=(0, 2))
    basis1 = pi_module.BasisAtom("Rb", n=(59, 63), l=(0, 1))

    eff_system = pi_module.EffectiveSystemPair([(ket, ket)])
    eff_system.set_basis_atoms(basis0)
    assert eff_system.basis_atoms[0] is basis0
    assert eff_system.basis_atoms[1] is basis0

    eff_system = pi_module.EffectiveSystemPair([(ket, ket)])
    eff_system.set_basis_atoms([basis0, basis1])
    assert eff_system.basis_atoms[0] is basis0
    assert eff_system.basis_atoms[1] is basis1

    eff_system = pi_module.EffectiveSystemPair([(ket, ket)])
    eff_system.basis_atoms  # noqa: B018
    # Now basis_atoms is set; trying to set it again should fail
    with pytest.raises(RuntimeError, match="basis_atoms"):
        eff_system.set_basis_atoms(basis0)


def test_set_system_atoms(pi_module: PairinteractionModule) -> None:
    ket = pi_module.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    basis0 = pi_module.BasisAtom("Rb", n=(58, 64), l=(0, 2))
    basis1 = pi_module.BasisAtom("Rb", n=(59, 63), l=(0, 1))
    system0 = pi_module.SystemAtom(basis0)
    system1 = pi_module.SystemAtom(basis1)

    eff_system = pi_module.EffectiveSystemPair([(ket, ket)])
    eff_system.set_system_atoms(system0)
    assert eff_system.system_atoms[0] is system0
    assert eff_system.system_atoms[1] is system0

    pi_module.diagonalize([system0, system1])
    eff_system = pi_module.EffectiveSystemPair([(ket, ket)])
    eff_system.set_system_atoms([system0, system1])
    assert eff_system.system_atoms[0] is system0
    assert eff_system.system_atoms[1] is system1

    eff_system = pi_module.EffectiveSystemPair([(ket, ket)])
    eff_system.basis_atoms  # noqa: B018
    with pytest.raises(RuntimeError, match="system_atoms"):
        eff_system.set_system_atoms(system0)

    eff_system = pi_module.EffectiveSystemPair([(ket, ket)])
    eff_system.set_electric_field([0, 0, 0], "volt/cm")
    with pytest.raises(RuntimeError, match="system atom parameters"):
        eff_system.set_system_atoms(system0)


def test_set_system_pair(pi_module: PairinteractionModule, system_pair_sample: SystemPair) -> None:
    ket = pi_module.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    basis = pi_module.BasisAtom("Rb", n=(58, 64), l=(0, 2))

    eff_system = pi_module.EffectiveSystemPair([(ket, ket)])
    eff_system.set_distance(10, unit="micrometer")
    with pytest.raises(RuntimeError, match="system pair parameters"):
        eff_system.set_system_pair(system_pair_sample)

    eff_system = pi_module.EffectiveSystemPair([(ket, ket)])
    eff_system.set_basis_atoms(basis)
    with pytest.raises(RuntimeError, match="electric_field"):
        eff_system.set_electric_field([0, 0, 1], "V/cm")


def test_set_distance_with_pint_quantity(pi_module: PairinteractionModule) -> None:
    """Test that set_distance accepts a pint.Quantity directly."""
    ket = pi_module.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    distance = 10 * ureg.micrometer
    eff_system = pi_module.EffectiveSystemPair([(ket, ket)])
    eff_system.set_distance(distance)
    assert eff_system._distance_vector is not None

    eff_system = pi_module.EffectiveSystemPair([(ket, ket)])
    eff_system.set_distance_vector([0, 0, distance])
    assert eff_system._distance_vector is not None


def test_create_basis_pair_with_delta_energy(pi_module: PairinteractionModule) -> None:
    """Test that create_basis_pair works with explicit delta_energy parameters."""
    ket1 = pi_module.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    ket2 = pi_module.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5)
    eff_system = pi_module.EffectiveSystemPair([(ket1, ket2)])
    eff_system.create_basis_pair(delta_energy=3.0, delta_energy_unit="GHz")
    assert eff_system._basis_pair is not None
    assert eff_system.basis_pair.number_of_kets > 0


def test_create_basis_pair_with_number_of_kets(pi_module: PairinteractionModule) -> None:
    """Test that create_basis_pair works with an explicit number_of_kets."""
    ket1 = pi_module.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    ket2 = pi_module.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5)
    eff_system = pi_module.EffectiveSystemPair([(ket1, ket2)])
    eff_system.create_basis_pair(number_of_kets=200)
    assert eff_system.basis_pair.number_of_kets >= 200


def test_get_pair_energies(pi_module: PairinteractionModule) -> None:
    """Test that get_pair_energies returns the correct number of energies."""
    ket1 = pi_module.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    ket2 = pi_module.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5)
    eff_system = pi_module.EffectiveSystemPair([(ket1, ket2), (ket2, ket1)])
    energies = eff_system.get_pair_energies(unit="GHz")
    assert len(energies) == 2
    # Both pair states (ket1,ket2) and (ket2,ket1) have the same energy
    assert np.isclose(energies[0], energies[1])
