# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix

    from pairinteraction.units import NDArray


logger = logging.getLogger(__name__)


def calculate_perturbative_hamiltonian(
    hamiltonian: csr_matrix,
    model_inds: list[int],
    perturbation_order: int,
) -> tuple[dict[int, NDArray], csr_matrix]:
    r"""Calculate the perturbative Hamiltonian up to a given order.

    This function calculates both the effective Hamiltonian, spanned up by the states of the model space,
    as well as the effective (perturbed) basisvectors due to interactions with the exterior space
    up to the desired order of perturbation theory.

    Args:
        hamiltonian: Hermitian matrix.
        model_inds: List of indices corresponding to the states that span up the model space.
        perturbation_order: Order up to which the perturbation theory is expanded.
            Support up to third order. Default is second order.

    Returns:
        Effective hamiltonians as dictionaries mapping perturbation orders to :math:`m \times m` matrices.
        Effective basisvectors in perturbation theory, returned as a sparse matrix. Each row represents one basisvector.

    """
    m_inds = np.asarray(model_inds, dtype=int)
    o_inds = np.setdiff1d(np.arange(hamiltonian.shape[0]), m_inds)
    eff_h_dict, eff_vecs = _calculate_unsorted_perturbative_hamiltonian(hamiltonian, m_inds, o_inds, perturbation_order)

    # resort eigvec to original order
    all_inds = np.append(m_inds, o_inds)
    all_inds_positions = np.argsort(all_inds)
    eff_vecs = eff_vecs[:, all_inds_positions]

    # include the hermitian conjugate part of the effective Hamiltonian
    for order, h_eff in eff_h_dict.items():
        eff_h_dict[order] = 0.5 * (h_eff + h_eff.conj().T)

    return eff_h_dict, eff_vecs


def _calculate_unsorted_perturbative_hamiltonian(
    hamiltonian: csr_matrix,
    m_inds: NDArray,
    o_inds: NDArray,
    perturbation_order: int,
) -> tuple[dict[int, NDArray], csr_matrix]:
    # This function is outsourced from calculate_perturbative_hamiltonian to allow for better type checking
    energies = np.real_if_close(hamiltonian.diagonal())
    if any(np.iscomplex(energies)):
        logger.error("The Hamiltonian has complex entries on the diagonal, this might lead to unexpected results.")

    energies_m = energies[m_inds]

    eff_h_dict: dict[int, NDArray] = {}  # perturbation_order -> h_eff

    eff_h_dict[0] = np.diag(energies_m)
    eff_vecs = sparse.csr_matrix(
        sparse.hstack(
            [sparse.eye(len(m_inds), len(m_inds), format="csr"), sparse.csr_matrix((len(m_inds), len(o_inds)))]
        )
    )

    if perturbation_order == 0:
        return eff_h_dict, eff_vecs

    v_offdiag = hamiltonian - sparse.diags(energies)
    v_mm = v_offdiag[np.ix_(m_inds, m_inds)]
    eff_h_dict[1] = v_mm.toarray()

    if perturbation_order == 1:
        return eff_h_dict, eff_vecs

    energies_diff = energies_m[np.newaxis, :] - energies[o_inds, np.newaxis]
    with np.errstate(divide="ignore"):
        delta_e_em = 1 / energies_diff
    v_me = v_offdiag[np.ix_(m_inds, o_inds)]
    eff_h_dict[2] = (v_me @ ((v_me.conj().T).multiply(delta_e_em))).toarray()

    addition_mm = sparse.csr_matrix((len(m_inds), len(m_inds)))
    addition_me = sparse.csr_matrix(((v_me.conj().T).multiply(delta_e_em)).T)
    if np.isinf(np.isinf(addition_me.data)).any():
        logger.critical(
            "Detected 'inf' entries in the effective basisvectors. "
            "This might happen, if you forgot to include a degenerate state in the model space. "
        )

    nan_idx = np.where(np.isnan(addition_me.data))[0]
    nonzero_inds = addition_me.nonzero()
    for i in nan_idx:
        value = addition_me.data[i]
        if np.isinf(value.real):
            row, col = nonzero_inds[0][i], nonzero_inds[1][i]
            addition_me[row, col] = value.real
        else:
            logger.error("Detected unexpected 'nan' entries in the effective basisvectors.")

    eff_vecs = eff_vecs + sparse.hstack([addition_mm, addition_me])

    if perturbation_order == 2:
        return eff_h_dict, eff_vecs

    diff = energies_m[np.newaxis, :] - energies_m[:, np.newaxis]
    diff = np.where(diff == 0, np.inf, diff)
    delta_e_mm = 1 / diff
    v_ee = v_offdiag[np.ix_(o_inds, o_inds)]
    if len(m_inds) > 1:
        logger.debug(
            "At third order, the effective basisvectors are currently only valid, "
            "when only one state is in the model space. "
            "Take care with interpretation of the effective basisvectors."
        )
    _intermediate_term = v_ee @ ((v_me.conj().T).multiply(delta_e_em)) - ((v_me.conj().T).multiply(delta_e_em)) @ v_mm
    eff_h_dict[3] = (v_me @ (_intermediate_term.multiply(delta_e_em))).toarray()

    addition_mm_diag = -0.5 * sparse.csr_matrix(
        sparse.diags((v_me @ ((v_me.conj().T).multiply(np.square(delta_e_em)))).diagonal())
    )
    addition_mm_offdiag = sparse.csr_matrix(((v_me @ (v_me.conj().T).multiply(delta_e_em)).multiply(delta_e_mm)).T)
    addition_me = sparse.csr_matrix(((v_ee @ ((v_me.conj().T).multiply(delta_e_em))).multiply(delta_e_em)).T)
    addition_me_2 = sparse.csr_matrix(((v_me.conj().T @ ((v_mm.conj().T).multiply(delta_e_mm))).multiply(delta_e_em)).T)
    eff_vecs = eff_vecs + sparse.hstack([addition_mm_diag + addition_mm_offdiag, addition_me + addition_me_2])

    if perturbation_order == 3:
        return eff_h_dict, eff_vecs

    raise ValueError("Perturbation theory is only implemented for orders [0, 1, 2, 3].")
