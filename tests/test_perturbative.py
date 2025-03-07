from typing import Union

import numpy as np
import pytest
from scipy import sparse

import pairinteraction.perturbative as perturbative
import pairinteraction.real as pi
from pairinteraction.perturbative.perturbative import _calculate_perturbative_hamiltonian


# define efficient equality function for sparse matrices
def _check_sparse_matrices_equal(A: sparse.csr_matrix, B: sparse.csr_matrix) -> Union[np.bool, bool]:
    """Check for equality of sparse matrices.

    This functions compares two sparse matrices in compressed sparse row format on their equality.

    Args:
        A: A sparse matrix in csr format.
        B: A sparse matrix in csr format.

    Returns:
        bool: True if matrices are equal, False if not.

    """
    A.sort_indices()
    B.sort_indices()
    if (
        A.format == "csr"
        and B.format == "csr"
        and len(A.indices) == len(B.indices)
        and len(A.indptr) == len(B.indptr)
        and len(A.data) == len(B.data)
    ):
        return (
            np.all(A.indices == B.indices)
            and np.all(A.indptr == B.indptr)
            and np.allclose(A.data, B.data, rtol=0, atol=1e-14)
        )
    else:
        return False


def _create_system_pair_sample():
    bases = [
        pi.BasisAtom(
            species="Rb",
            n=(59, 63),
            l=(0, 1),
            m=(-1.5, 1.5),
        )
        for i in [1, 2]
    ]
    systems = [pi.SystemAtom(basis=basis) for basis in bases]
    for system in systems:
        system.enable_diamagnetism(False)
        system.set_magnetic_field([0, 0, 1e-3], "gauss")
    if not all(system.is_diagonal for system in systems):
        pi.diagonalize(systems, diagonalizer="eigen", sort_by_energy=False)
    basis_pair = pi.BasisPair(systems)
    system_pair = pi.SystemPair(basis_pair)
    theta = 0
    r = 12
    system_pair.set_distance_vector(r * np.array([np.sin(theta), 0, np.cos(theta)]), "micrometer")
    system_pair.set_order(3)
    return system_pair


# test for mathematic functionality
def test_perturbative_calculation():
    H = sparse.csr_matrix([[0, 1, 1, 0, 2], [1, 1, 0, 1, 3], [1, 0, 10, 1, 0], [0, 1, 1, 11, 1], [2, 3, 0, 1, 12]])
    model_space_indices = np.array([0, 1])

    H_eff, eig_perturb = _calculate_perturbative_hamiltonian(H=H, model_space_indices=model_space_indices, order=0)

    assert np.any(H_eff == np.array([[0, 0], [0, 1]]))
    assert _check_sparse_matrices_equal(eig_perturb, sparse.eye(2, 5, k=0, format="csr"))

    H_eff, eig_perturb = _calculate_perturbative_hamiltonian(H=H, model_space_indices=model_space_indices, order=1)

    assert np.any(H_eff == np.array([[0, 1], [1, 1]]))
    assert _check_sparse_matrices_equal(eig_perturb, sparse.csr_matrix([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]]))

    H_eff, eig_perturb = _calculate_perturbative_hamiltonian(H=H, model_space_indices=model_space_indices, order=2)
    a_00 = 0 + 1 * 1 / (0 - 10) + 2 * 2 / (0 - 12)
    a_01 = 1 + 2 * 3 / (0 - 12)
    a_10 = 1 + 3 * 2 / (1 - 12)
    a_11 = 1 + 1 * 1 / (1 - 11) + 3 * 3 / (1 - 12)
    H_new = np.array([[a_00, (a_01 + a_10) / 2], [(a_01 + a_10) / 2, a_11]])
    assert np.any(H_eff == H_new)

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

    H_eff, eig_perturb = _calculate_perturbative_hamiltonian(H=H, model_space_indices=model_space_indices, order=3)
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
    H_new = sparse.csr_matrix([[a_00, (a_01 + a_10) / 2], [(a_01 + a_10) / 2, a_11]])
    assert np.any(H_eff == H_new)

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


# check that dispersion coefficients are correctly calculated


def test_c3_coefficient():
    system_pair = _create_system_pair_sample()

    ket_tuple_list = [
        (pi.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5), pi.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5)),
        (pi.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5), pi.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)),
    ]
    C3 = perturbative.get_c3_from_system(
        system_pair=system_pair, ket_tuple_list=ket_tuple_list, unit="planck_constant * gigahertz * micrometer^3"
    )
    assert np.isclose(-0.5 * C3.magnitude, 3.1515)


# # calculate and check that C6 coefficient is correct


def test_c6_coefficient():
    system_pair = _create_system_pair_sample()
    ket_atom = pi.KetAtom(species="Rb", n=61, l=0, j=0.5, m=0.5)
    C6 = perturbative.get_c6_from_system(
        ket_tuple=(ket_atom, ket_atom), system_pair=system_pair, unit="planck_constant * gigahertz * micrometer^6"
    )
    assert np.isclose(C6.magnitude, 167.92)


# check that resonance is correctly detected


def test_resonance_detection():
    system_pair = _create_system_pair_sample()
    ket_tuple_list = [
        (pi.KetAtom(species="Rb", n=61, l=0, j=0.5, m=0.5), pi.KetAtom(species="Rb", n=61, l=1, j=0.5, m=0.5))
    ]
    with pytest.warns(RuntimeWarning, match="divide by zero encountered in divide"):
        with pytest.raises(ValueError):
            H_eff, _ = perturbative.get_effective_hamiltonian_from_system(
                ket_tuple_list=ket_tuple_list, system_pair=system_pair, order=2, required_overlap=0.8
            )
