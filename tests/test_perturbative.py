# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging

import numpy as np
import pairinteraction.real as pi
import pytest
from pairinteraction import perturbative
from pairinteraction.perturbative.perturbative import _calculate_perturbative_hamiltonian
from scipy import sparse


# ---------------------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------------------
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


def _create_system_pair_sample() -> pi.SystemPair:
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


# ---------------------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------------------
def test_perturbative_calculation1(caplog: pytest.LogCaptureFixture) -> None:
    """Test of mathematic functionality."""
    hamiltonian = sparse.csr_matrix(
        [[0, 1, 1, 0, 2], [1, 1, 0, 1, 3], [1, 0, 10, 1, 0], [0, 1, 1, 11, 1], [2, 3, 0, 1, 12]]
    )
    model_space_indices = [0, 1]

    # Order 0
    hamiltonian_eff, eig_perturb = _calculate_perturbative_hamiltonian(hamiltonian, model_space_indices, order=0)
    assert np.any(hamiltonian_eff == np.array([[0, 0], [0, 1]]))
    assert _check_sparse_matrices_equal(eig_perturb, sparse.csr_matrix(sparse.eye(2, 5, k=0)))

    # Order 1
    hamiltonian_eff, eig_perturb = _calculate_perturbative_hamiltonian(hamiltonian, model_space_indices, order=1)
    assert np.any(hamiltonian_eff == np.array([[0, 1], [1, 1]]))
    assert _check_sparse_matrices_equal(eig_perturb, sparse.csr_matrix([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]]))

    # Order 2
    hamiltonian_eff, eig_perturb = _calculate_perturbative_hamiltonian(hamiltonian, model_space_indices, order=2)
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
    assert _check_sparse_matrices_equal(eig_perturb, sparse.csr_matrix(sparse.vstack([v0, v1])))

    # Order 3
    with caplog.at_level(logging.ERROR):
        hamiltonian_eff, eig_perturb = _calculate_perturbative_hamiltonian(hamiltonian, model_space_indices, order=3)
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
    assert _check_sparse_matrices_equal(eig_perturb, sparse.csr_matrix(sparse.vstack([v0, v1])))


def test_c3_coefficient() -> None:
    """Test whether dispersion coefficients are correctly calculated."""
    system_pair = _create_system_pair_sample()

    ket_tuple_list = [
        (pi.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5), pi.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5)),
        (pi.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5), pi.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)),
    ]
    c3 = perturbative.get_c3_from_system(
        system_pair=system_pair, ket_tuple_list=ket_tuple_list, unit="planck_constant * gigahertz * micrometer^3"
    )
    assert np.isclose(-0.5 * c3, 3.1515)


def test_c6_coefficient() -> None:
    """Test whether the c6 coefficient is correct."""
    system_pair = _create_system_pair_sample()
    ket_atom = pi.KetAtom(species="Rb", n=61, l=0, j=0.5, m=0.5)
    c6 = perturbative.get_c6_from_system(
        ket_tuple=(ket_atom, ket_atom), system_pair=system_pair, unit="planck_constant * gigahertz * micrometer^6"
    )
    assert np.isclose(c6, 167.92)


def test_resonance_detection(caplog: pytest.LogCaptureFixture) -> None:
    """Test whether resonance is correctly detected."""
    system_pair = _create_system_pair_sample()
    ket_tuple_list = [
        (pi.KetAtom(species="Rb", n=61, l=0, j=0.5, m=0.5), pi.KetAtom(species="Rb", n=61, l=1, j=0.5, m=0.5))
    ]
    with (
        pytest.warns(RuntimeWarning, match="divide by zero encountered in divide"),
        pytest.raises(ValueError, match="Perturbative Calculation not possible due to resonances."),
        caplog.at_level(logging.CRITICAL),
    ):
        perturbative.get_effective_hamiltonian_from_system(
            ket_tuple_list=ket_tuple_list, system_pair=system_pair, order=2, required_overlap=0.8
        )
