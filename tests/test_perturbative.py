# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging

import numpy as np
import pairinteraction.real as pi
import pytest
from pairinteraction import perturbative
from pairinteraction.perturbative.effective_hamiltonian import _calculate_perturbative_hamiltonian
from pairinteraction.units import ureg
from scipy import sparse


def _check_sparse_matrices_equal(matrix_a: sparse.csr_matrix, matrix_b: sparse.csr_matrix) -> bool:
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
    if (
        matrix_a.format == "csr"
        and matrix_b.format == "csr"
        and len(matrix_a.indices) == len(matrix_b.indices)
        and len(matrix_a.indptr) == len(matrix_b.indptr)
        and len(matrix_a.data) == len(matrix_b.data)
    ):
        return bool(
            np.all(matrix_a.indices == matrix_b.indices)
            and np.all(matrix_a.indptr == matrix_b.indptr)
            and np.allclose(matrix_a.data, matrix_b.data, rtol=0, atol=1e-14)
        )
    return False


@pytest.fixture
def system_pair() -> pi.SystemPair:
    basis = pi.BasisAtom(
        species="Rb",
        n=(59, 63),
        l=(0, 1),
        m=(-1.5, 1.5),
    )
    system = pi.SystemAtom(basis=basis)
    system.set_diamagnetism_enabled(False)
    system.set_magnetic_field([0, 0, 1e-3], "gauss")
    if not system.is_diagonal:
        pi.diagonalize([system], diagonalizer="eigen", sort_by_energy=False)
    basis_pair = pi.BasisPair([system, system])
    system_pair = pi.SystemPair(basis_pair)
    theta = 0
    r = 12
    system_pair.set_distance_vector(r * np.array([np.sin(theta), 0, np.cos(theta)]), "micrometer")
    system_pair.set_interaction_order(3)
    return system_pair


def test_perturbative_calculation(caplog) -> None:
    """Test of mathematic functionality."""
    H = sparse.csr_matrix([[0, 1, 1, 0, 2], [1, 1, 0, 1, 3], [1, 0, 10, 1, 0], [0, 1, 1, 11, 1], [2, 3, 0, 1, 12]])
    model_space_indices = [0, 1]

    hamiltonian_eff, eig_perturb = _calculate_perturbative_hamiltonian(H, model_space_indices, perturbation_order=0)

    assert np.any(hamiltonian_eff == np.array([[0, 0], [0, 1]]))
    assert _check_sparse_matrices_equal(eig_perturb, sparse.eye(2, 5, k=0, format="csr"))

    hamiltonian_eff, eig_perturb = _calculate_perturbative_hamiltonian(H, model_space_indices, perturbation_order=1)

    assert np.any(hamiltonian_eff == np.array([[0, 1], [1, 1]]))
    assert _check_sparse_matrices_equal(eig_perturb, sparse.csr_matrix([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]]))

    hamiltonian_eff, eig_perturb = _calculate_perturbative_hamiltonian(H, model_space_indices, perturbation_order=2)
    a_00 = 0 + 1 * 1 / (0 - 10) + 2 * 2 / (0 - 12)
    a_01 = 1 + 2 * 3 / (0 - 12)
    a_10 = 1 + 3 * 2 / (1 - 12)
    a_11 = 1 + 1 * 1 / (1 - 11) + 3 * 3 / (1 - 12)
    hamiltonian_new = np.array([[a_00, (a_01 + a_10) / 2], [(a_01 + a_10) / 2, a_11]])
    assert np.any(hamiltonian_eff == hamiltonian_new)

    v0 = (
        sparse.eye(1, 5, k=0, format="csr")
        + 1 / (0 - 10) * sparse.eye(1, 5, k=2, format="csr")
        + 2 / (0 - 12) * sparse.eye(1, 5, k=4, format="csr")
    )
    v1 = (
        sparse.eye(1, 5, k=1, format="csr")
        + 1 / (1 - 11) * sparse.eye(1, 5, k=3, format="csr")
        + 3 / (1 - 12) * sparse.eye(1, 5, k=4, format="csr")
    )
    assert _check_sparse_matrices_equal(eig_perturb, sparse.vstack([v0, v1]))

    with caplog.at_level(logging.ERROR):
        hamiltonian_eff, eig_perturb = _calculate_perturbative_hamiltonian(H, model_space_indices, perturbation_order=3)
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
    hamiltonian_new = sparse.csr_matrix([[a_00, (a_01 + a_10) / 2], [(a_01 + a_10) / 2, a_11]])
    assert np.any(hamiltonian_eff == hamiltonian_new)

    v0 += (
        -0.5 * (1 * 1 / (0 - 10) ** 2 + 2 * 2 / (0 - 12) ** 2) * sparse.eye(1, 5, k=0, format="csr")
        + (1 * 1 / ((0 - 10) * (0 - 11)) + 2 * 1 / ((0 - 11) * (0 - 12))) * sparse.eye(1, 5, k=3, format="csr")
        + 1 * 1 / ((0 - 11) * (0 - 1)) * sparse.eye(1, 5, k=3, format="csr")
        + 3 * 1 / ((0 - 1) * (0 - 12)) * sparse.eye(1, 5, k=4, format="csr")
        + 3 * 2 / ((0 - 1) * (0 - 12)) * sparse.eye(1, 5, k=1, format="csr")
    )
    v1 += (
        -0.5 * (1 * 1 / (1 - 11) ** 2 + 3 * 3 / (1 - 12) ** 2) * sparse.eye(1, 5, k=1, format="csr")
        + 1 * 1 / ((1 - 10) * (1 - 11)) * sparse.eye(1, 5, k=2, format="csr")
        + 1 * 3 / ((1 - 11) * (1 - 12)) * sparse.eye(1, 5, k=3, format="csr")
        + 1 * 1 / ((1 - 11) * (1 - 12)) * sparse.eye(1, 5, k=4, format="csr")
        + 1 * 1 / ((1 - 0) * (1 - 10)) * sparse.eye(1, 5, k=2, format="csr")
        + 2 * 1 / ((1 - 0) * (1 - 12)) * sparse.eye(1, 5, k=4, format="csr")
        + 3 * 2 / ((1 - 0) * (1 - 12)) * sparse.eye(1, 5, k=0, format="csr")
    )
    assert _check_sparse_matrices_equal(eig_perturb, sparse.vstack([v0, v1]))


def test_c3_with_system(system_pair: pi.SystemPair) -> None:
    """Test whether the C3 coefficient with a given system is calculated correctly."""
    ket_tuple_list = [
        (pi.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5), pi.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5)),
        (pi.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5), pi.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)),
    ]
    c3 = perturbative.get_c3_from_system(ket_tuple_list, system_pair, unit="planck_constant * gigahertz * micrometer^3")
    assert np.isclose(-0.5 * c3, 3.1515)


def test_c3_without_system() -> None:
    """Test whether the C3 coefficient with automatically constructed system is calculated correctly."""
    ket_tuple_list = [
        (pi.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5), pi.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5)),
        (pi.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5), pi.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)),
    ]
    magnetic_field = ureg.Quantity([0, 0, 10], "gauss")
    electric_field = ureg.Quantity([0, 0, 0], "volt/cm")
    distance_vector = ureg.Quantity([0, 0, 500], "micrometer")
    c3 = perturbative.get_c3(
        ket_tuple_list,
        distance_vector,
        electric_field,
        magnetic_field,
        unit="planck_constant * gigahertz * micrometer^3",
    )
    assert np.isclose(-0.5 * c3, 3.2188)


def test_c6_with_system(system_pair: pi.SystemPair) -> None:
    """Test whether the C6 coefficient with a given system is calculated correctly."""
    ket_atom = pi.KetAtom(species="Rb", n=61, l=0, j=0.5, m=0.5)
    c6 = perturbative.get_c6_from_system(
        (ket_atom, ket_atom), system_pair, unit="planck_constant * gigahertz * micrometer^6"
    )
    assert np.isclose(c6, 167.92)


def test_c6_without_system() -> None:
    """Test whether the C6 coefficient with automatically constructed system is calculated correctly."""
    magnetic_field = ureg.Quantity([0, 0, 10], "gauss")
    electric_field = ureg.Quantity([0, 0, 0], "volt/cm")
    distance_vector = ureg.Quantity([0, 0, 500], "micrometer")
    ket_atom = pi.KetAtom(species="Rb", n=61, l=0, j=0.5, m=0.5)
    c6 = perturbative.get_c6(
        (ket_atom, ket_atom),
        distance_vector,
        electric_field,
        magnetic_field,
        unit="planck_constant * gigahertz * micrometer^6",
    )
    assert np.isclose(c6, 169.189)


def test_resonance_detection(system_pair: pi.SystemPair, caplog) -> None:
    """Test whether resonance is correctly detected."""
    ket_tuple_list = [
        (pi.KetAtom(species="Rb", n=61, l=0, j=0.5, m=0.5), pi.KetAtom(species="Rb", n=61, l=1, j=0.5, m=0.5))
    ]
    with (
        caplog.at_level(logging.CRITICAL),
        pytest.warns(RuntimeWarning, match="invalid value encountered in multiply"),
    ):
        perturbative.get_effective_hamiltonian_from_system(ket_tuple_list, system_pair, perturbation_order=2)
