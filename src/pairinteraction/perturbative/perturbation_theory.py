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


def calculate_perturbative_hamiltonian(
    hamiltonian: "csr_matrix",
    model_inds: list[int],
    perturbation_order: int,
) -> tuple[dict[int, "NDArray"], "csr_matrix"]:
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
    eff_h_dict, eff_eigvec = _calculate_unsorted_perturbative_hamiltonian(
        hamiltonian, m_inds, o_inds, perturbation_order
    )

    # resort eigvec to original order
    all_inds = np.append(m_inds, o_inds)
    all_inds_positions = np.argsort(all_inds)
    eff_eigvec = eff_eigvec[:, all_inds_positions]

    # include the hermitian conjugate part of the effective Hamiltonian
    for order, h_eff in eff_h_dict.items():
        eff_h_dict[order] = 0.5 * (h_eff + h_eff.conj().T)

    return eff_h_dict, eff_eigvec


def _calculate_unsorted_perturbative_hamiltonian(
    hamiltonian: "csr_matrix",
    m_inds: "NDArray",
    o_inds: "NDArray",
    perturbation_order: int,
) -> tuple[dict[int, "NDArray"], "csr_matrix"]:
    # This function is outsourced from calculate_perturbative_hamiltonian to allow for better type checking
    h0 = hamiltonian.diagonal()
    h0_m = h0[m_inds]

    eff_h_dict: dict[int, NDArray] = {}  # perturbation_order -> h_eff

    eff_h_dict[0] = np.diag(h0_m)
    eff_eigvec = sparse.hstack(
        [sparse.eye(len(m_inds), len(m_inds), format="csr").tocsr(), sparse.csr_matrix((len(m_inds), len(o_inds)))]
    ).tocsr()

    if perturbation_order == 0:
        return eff_h_dict, eff_eigvec

    v_offdiag = hamiltonian - sparse.diags(h0)
    v_mm = v_offdiag[np.ix_(m_inds, m_inds)]
    eff_h_dict[1] = v_mm.toarray()

    if perturbation_order == 1:
        return eff_h_dict, eff_eigvec

    h0_e = h0[o_inds]
    v_me = v_offdiag[np.ix_(m_inds, o_inds)]
    with np.errstate(divide="ignore"):
        delta_e_em = 1 / (h0_m[np.newaxis, :] - h0_e[:, np.newaxis])
    eff_h_dict[2] = (v_me @ ((v_me.conj().T).multiply(delta_e_em))).toarray()

    addition_mm = sparse.csr_matrix((len(m_inds), len(m_inds)))
    addition_me = sparse.csr_matrix(((v_me.conj().T).multiply(delta_e_em)).T)
    if np.isinf(addition_me.data).any():
        logger.critical(
            "Detected 'inf' entries in the effective basisvectors. "
            "This might happen, if you forgot to include a degenerate state in the model space. "
        )
    eff_eigvec = eff_eigvec + sparse.hstack([addition_mm, addition_me])

    if perturbation_order == 2:
        return eff_h_dict, eff_eigvec

    diff = h0_m[np.newaxis, :] - h0_m[:, np.newaxis]
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
    eff_eigvec = eff_eigvec + sparse.hstack([addition_mm_diag + addition_mm_offdiag, addition_me + addition_me_2])

    if perturbation_order == 3:
        return eff_h_dict, eff_eigvec

    raise ValueError("Perturbation theory is only implemented for orders [0, 1, 2, 3].")
