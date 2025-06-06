# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from typing import TYPE_CHECKING

import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix

    from pairinteraction.units import NDArray


logger = logging.getLogger(__name__)


def _calculate_perturbative_hamiltonian(
    hamiltonian: "csr_matrix", model_inds: list[int], order: int = 2, return_only_specified_order: bool = False
) -> tuple["NDArray", "csr_matrix"]:
    r"""Calculate the perturbative Hamiltonian at a given order.

    This function takes a Hamiltonian as a sparse matrix which is diagonal in the unperturbed basis
    and list of indices spanning up the model space.
    It calculates both the effective Hamiltonian, spanned up by the states of the model space, as well as the
    perturbed eigenstates due to interactions with the exterior space in the desired order of perturbation theory.
    The output is either the full Hamiltonian up to the order of perturbation theory, or only the corrections at
    a given order.

    Args:
        hamiltonian: Quadratic hermitian matrix. Perturbative terms are assumed to be only off-diagonal.
        model_inds: List of indices corresponding to the states that span up the model space.
        order: Order up to which the perturbation theory is expanded. Support up to third order.
            Default is second order.
        return_only_specified_order: If True, the returned effective Hamiltonian will only contain the specified order.
            Default is False, which returns the sum of all orders up to the specified order.

    Returns:
        Effective Hamiltonian as a :math:`m \times m` matrix, where m is the length of `ket_tuple_list`
        Eigenvectors in perturbation theory due to interaction with states out of the model
            space, returned as a sparse matrix in compressed row format. Each row represent the
            corresponding eigenvector

    """
    m_inds = np.array(model_inds)
    o_inds = np.setdiff1d(np.arange(hamiltonian.shape[0]), m_inds)
    h_eff, eigvec_perturb = _calculate_unsorted_perturbative_hamiltonian(
        hamiltonian, m_inds, o_inds, order, return_only_specified_order
    )

    # resort eigvec to original order
    all_inds = np.append(m_inds, o_inds)
    all_inds_positions = np.argsort(all_inds)
    eigvec_perturb = eigvec_perturb[:, all_inds_positions]

    # include the hermitian conjugate part of the effective Hamiltonian
    h_eff = 0.5 * (h_eff + h_eff.conj().T)

    return h_eff, eigvec_perturb


def _calculate_unsorted_perturbative_hamiltonian(  # noqa: C901
    hamiltonian: "csr_matrix",
    m_inds: "NDArray",
    o_inds: "NDArray",
    order: int,
    return_only_specified_order: bool = False,
) -> tuple["NDArray", "csr_matrix"]:
    # This function is outsourced from _calculate_perturbative_hamiltonian to allow for better type checking
    if order not in [0, 1, 2, 3]:
        raise ValueError("Perturbation theory is only implemented for orders [0, 1, 2, 3].")

    h0 = hamiltonian.diagonal()
    h0_m = h0[m_inds]
    h_eff = np.zeros_like(np.diag(h0_m))

    if not return_only_specified_order or order == 0:
        h_eff += np.diag(h0_m)
    eigvec_perturb = sparse.hstack(
        [sparse.eye(len(m_inds), len(m_inds), format="csr").tocsr(), sparse.csr_matrix((len(m_inds), len(o_inds)))]
    ).tocsr()

    if order < 1:
        return h_eff, eigvec_perturb

    v_offdiag = hamiltonian - sparse.diags(h0)
    v_mm = v_offdiag[np.ix_(m_inds, m_inds)]
    if not return_only_specified_order or order == 1:
        h_eff += v_mm

    if order < 2:
        return h_eff, eigvec_perturb

    h0_e = h0[o_inds]
    v_me = v_offdiag[np.ix_(m_inds, o_inds)]
    with np.errstate(divide="ignore"):
        delta_e_em = 1 / (h0_m[np.newaxis, :] - h0_e[:, np.newaxis])
    if not return_only_specified_order or order == 2:
        h_eff += v_me @ ((v_me.conj().T).multiply(delta_e_em))
    addition_mm = sparse.csr_matrix((len(m_inds), len(m_inds)))
    addition_me = sparse.csr_matrix(((v_me.conj().T).multiply(delta_e_em)).T)
    eigvec_perturb = eigvec_perturb + sparse.hstack([addition_mm, addition_me])

    if order < 3:
        return h_eff, eigvec_perturb

    diff = h0_m[np.newaxis, :] - h0_m[:, np.newaxis]
    diff = np.where(diff == 0, np.inf, diff)
    delta_e_mm = 1 / diff
    v_ee = v_offdiag[np.ix_(o_inds, o_inds)]
    if len(m_inds) > 1:
        logger.warning(
            "At third order, the eigenstates are currently only valid when only one state is in the model space. "
            "Take care with interpreation of the perturbed eigenvectors."
        )
    if not return_only_specified_order or order == 3:
        h_eff += v_me @ (
            (v_ee @ ((v_me.conj().T).multiply(delta_e_em)) - ((v_me.conj().T).multiply(delta_e_em)) @ v_mm).multiply(
                delta_e_em
            )
        )
    addition_mm_diag = -0.5 * sparse.csr_matrix(
        sparse.diags((v_me @ ((v_me.conj().T).multiply(np.square(delta_e_em)))).diagonal())
    )
    addition_mm_offdiag = sparse.csr_matrix(((v_me @ (v_me.conj().T).multiply(delta_e_em)).multiply(delta_e_mm)).T)
    addition_me = sparse.csr_matrix(((v_ee @ ((v_me.conj().T).multiply(delta_e_em))).multiply(delta_e_em)).T)
    addition_me_2 = sparse.csr_matrix(((v_me.conj().T @ ((v_mm.conj().T).multiply(delta_e_mm))).multiply(delta_e_em)).T)
    eigvec_perturb = eigvec_perturb + sparse.hstack(
        [addition_mm_diag + addition_mm_offdiag, addition_me + addition_me_2]
    )

    if order < 4:
        return h_eff, eigvec_perturb

    raise ValueError("Perturbation theory is only implemented for orders [0, 1, 2, 3].")
