# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from collections.abc import Collection
from typing import TYPE_CHECKING, Any, Optional, Union, overload

import numpy as np
from scipy import sparse

from pairinteraction.perturbative.create_system import create_system_for_perturbative
from pairinteraction.units import QuantityArray

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix

    from pairinteraction._wrapped.basis.basis_pair import BasisPair
    from pairinteraction._wrapped.ket.ket_pair import KetAtomTuple
    from pairinteraction._wrapped.system.system_pair import SystemPair
    from pairinteraction.units import NDArray, PintArray, PintArrayLike


logger = logging.getLogger(__name__)


@overload
def get_effective_hamiltonian_from_system(
    ket_tuples: Collection["KetAtomTuple"],
    system_pair: "SystemPair[Any]",
    perturbation_order: int = 2,
    *,
    return_only_specified_order: bool = False,
    unit: None = None,
) -> tuple["PintArray", "csr_matrix"]: ...


@overload
def get_effective_hamiltonian_from_system(
    ket_tuples: Collection["KetAtomTuple"],
    system_pair: "SystemPair[Any]",
    perturbation_order: int = 2,
    *,
    return_only_specified_order: bool = False,
    unit: str,
) -> tuple["NDArray", "csr_matrix"]: ...


def get_effective_hamiltonian_from_system(
    ket_tuples: Collection["KetAtomTuple"],
    system_pair: "SystemPair[Any]",
    perturbation_order: int = 2,
    return_only_specified_order: bool = False,
    unit: Optional[str] = None,
) -> tuple[Union["NDArray", "PintArray"], "csr_matrix"]:
    r"""Get the effective Hamiltonian at the desired order in Rayleigh-Schrödinger perturbation theory.

    This function takes a list of KetAtom 2-tuples (i.e. [(ketA_1, ketB_1), (ketA_2, ketB_2), ...]),
    which form the basis of the model space.
    The effective Hamiltonian is then calculated for this model space.
    The whole Hamiltonian is taken from a SystemPair object.
    Also check for resonances between states in the model space and other dipole coupled states.

    Args:
        ket_tuples: List of KetAtom 2-tuples that form the model space.
            The effective Hamiltonian is calculated for this model space.
        system_pair: SystemPair object that defines the Hamiltonian
            (including electric and magnetic fields, the interatomic distance, ...)
        perturbation_order: The order up to which the perturbation theory is taken into account.
            Support up to third order. Default is second order.
        return_only_specified_order: If True, the returned effective Hamiltonian will only contain the specified order.
            Default is False, which returns the sum of all orders up to the specified order.
        unit: The unit in which the effective hamiltonian will be returned.
            Default None will return a pint quantity.

    Returns:
        - Effective Hamiltonian as a :math:`m \times m` matrix, where m is the length of `ket_tuples`.
        - Perturbed Eigenvectors, returned as a sparse matrix.
          Each row corresponds to one eigenvector in the basis of the SystemPair.

    """
    if np.isinf(system_pair.get_distance().magnitude):
        raise ValueError(
            "The effective Hamiltonian can only be calculated for a given finite distance. "
            "If you want to calculate the C3 or C6 coefficient you still have to set some distance for the system "
            "or you can use the `get_c3(...)` or `get_c6(...)` functions."
        )

    model_inds = _get_model_inds(ket_tuples, system_pair)
    h_au = system_pair.get_hamiltonian().to_base_units().magnitude  # Hamiltonian in atomic units
    h_eff_au, eigvec_perturb = _calculate_perturbative_hamiltonian(
        h_au, model_inds, perturbation_order, return_only_specified_order=return_only_specified_order
    )

    check_for_resonances(ket_tuples, system_pair, eigvec_perturb, required_overlap=0.9)

    h_eff = QuantityArray.from_base_unit(h_eff_au, "energy").to_pint_or_unit(unit)
    return h_eff, eigvec_perturb


@overload
def get_effective_hamiltonian(
    ket_tuples: Collection["KetAtomTuple"],
    distance_vector: "PintArrayLike",
    electric_field: Optional["PintArrayLike"] = None,
    magnetic_field: Optional["PintArrayLike"] = None,
    multipole_order: int = 3,
    with_diamagnetism: bool = False,
    perturbation_order: int = 2,
    *,
    unit: None = None,
) -> tuple["PintArray", sparse.csr_matrix]: ...


@overload
def get_effective_hamiltonian(
    ket_tuples: Collection["KetAtomTuple"],
    distance_vector: "PintArrayLike",
    electric_field: Optional["PintArrayLike"] = None,
    magnetic_field: Optional["PintArrayLike"] = None,
    multipole_order: int = 3,
    with_diamagnetism: bool = False,
    perturbation_order: int = 2,
    *,
    unit: str,
) -> tuple["NDArray", sparse.csr_matrix]: ...


def get_effective_hamiltonian(
    ket_tuples: Collection["KetAtomTuple"],
    distance_vector: "PintArrayLike",
    electric_field: Optional["PintArrayLike"] = None,
    magnetic_field: Optional["PintArrayLike"] = None,
    multipole_order: int = 3,
    with_diamagnetism: bool = False,
    perturbation_order: int = 2,
    unit: Optional[str] = None,
) -> tuple[Union["NDArray", "PintArray"], sparse.csr_matrix]:
    r"""Get the perturbative Hamiltonian at a desired order in Rayleigh-Schrödinger perturbation theory.

    This function takes a list of KetAtom pairs (stored as tuple),
    which forms the basis of the model space for which the effective Hamiltonian is calculated.
    In addition it also checks for resonances between states in the model space and other dipole coupled states.

    Args:
        ket_tuples: List of tuples of KetAtom pairs, which form the basis of the model space.
        distance_vector: distance vector between the atoms.
        electric_field: electric field in the system.
        magnetic_field: magnetic field in the system.
        multipole_order: multipole-order of the interaction. Default is 3 (dipole-dipole).
        with_diamagnetism: True if diamagnetic term should be considered. Default is False.
        perturbation_order: order of perturbative calculation the system shall be used for. Default is 2.
        unit: The unit in which the effective Hamiltonian should be returned.
            Default None will return a pint quantity.

    Returns:
        - Effective Hamiltonian as a :math:`m \times m` matrix, where m is the length of `ket_tuples`,
          either in the given unit or as pint quantity if unit is None.
        - Eigenvectors in perturbation theory due to interaction with states out of the model space, returned as
          a sparse matrix in compressed row format. Each row represents the corresponding eigenvector.

    """
    system_pair = create_system_for_perturbative(
        ket_tuples,
        distance_vector,
        electric_field,
        magnetic_field,
        multipole_order,
        with_diamagnetism,
        perturbation_order,
    )
    return get_effective_hamiltonian_from_system(
        ket_tuples,
        system_pair,
        perturbation_order,
        unit=unit,
    )


def _calculate_perturbative_hamiltonian(
    hamiltonian: "csr_matrix",
    model_inds: list[int],
    perturbation_order: int = 2,
    *,
    return_only_specified_order: bool = False,
) -> tuple["NDArray", "csr_matrix"]:
    """Calculate the perturbative Hamiltonian at a given order."""
    m_inds = np.array(model_inds)
    o_inds = np.setdiff1d(np.arange(hamiltonian.shape[0]), m_inds)
    h_eff, eigvec_perturb = _calculate_unsorted_perturbative_hamiltonian(
        hamiltonian, m_inds, o_inds, perturbation_order, return_only_specified_order
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
    perturbation_order: int,
    return_only_specified_order: bool = False,
) -> tuple["NDArray", "csr_matrix"]:
    # This function is outsourced from _calculate_perturbative_hamiltonian to allow for better type checking
    if perturbation_order not in [0, 1, 2, 3]:
        raise ValueError("Perturbation theory is only implemented for orders [0, 1, 2, 3].")

    h0 = hamiltonian.diagonal()
    h0_m = h0[m_inds]
    h_eff = np.zeros_like(np.diag(h0_m))

    if not return_only_specified_order or perturbation_order == 0:
        h_eff += np.diag(h0_m)
    eigvec_perturb = sparse.hstack(
        [sparse.eye(len(m_inds), len(m_inds), format="csr").tocsr(), sparse.csr_matrix((len(m_inds), len(o_inds)))]
    ).tocsr()

    if perturbation_order < 1:
        return h_eff, eigvec_perturb

    v_offdiag = hamiltonian - sparse.diags(h0)
    v_mm = v_offdiag[np.ix_(m_inds, m_inds)]
    if not return_only_specified_order or perturbation_order == 1:
        h_eff += v_mm

    if perturbation_order < 2:
        return h_eff, eigvec_perturb

    h0_e = h0[o_inds]
    v_me = v_offdiag[np.ix_(m_inds, o_inds)]
    with np.errstate(divide="ignore"):
        delta_e_em = 1 / (h0_m[np.newaxis, :] - h0_e[:, np.newaxis])
    if not return_only_specified_order or perturbation_order == 2:
        h_eff += v_me @ ((v_me.conj().T).multiply(delta_e_em))
    addition_mm = sparse.csr_matrix((len(m_inds), len(m_inds)))
    addition_me = sparse.csr_matrix(((v_me.conj().T).multiply(delta_e_em)).T)
    eigvec_perturb = eigvec_perturb + sparse.hstack([addition_mm, addition_me])

    if perturbation_order < 3:
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
    if not return_only_specified_order or perturbation_order == 3:
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

    if perturbation_order < 4:
        return h_eff, eigvec_perturb

    raise ValueError("Perturbation theory is only implemented for orders [0, 1, 2, 3].")


def _get_model_inds(
    ket_tuples: Collection["KetAtomTuple"], system_pair: "SystemPair[BasisPair[Any, Any]]"
) -> list[int]:
    """Get the indices of all ket tuples in the basis of the system_pair."""
    model_inds = []
    for kets in ket_tuples:
        overlap = system_pair.basis.get_overlaps(kets)
        index = np.argmax(overlap)
        if overlap[index] == 0:
            raise ValueError(f"The pairstate {kets} is not part of the basis of the pair system.")
        if overlap[index] < 0.5:
            raise ValueError(f"The pairstate {kets} cannot be identified uniquely (max overlap: {overlap[index]}).")
        model_inds.append(int(index))
    return model_inds


def check_for_resonances(
    ket_tuples: Collection["KetAtomTuple"],
    system_pair: "SystemPair[BasisPair[Any, Any]]",
    eigvec_perturb: "csr_matrix",
    required_overlap: float,
) -> None:
    r"""Check for resonance between the states in the model space and other states.

    This function takes the perturbed eigenvectors of the perturbation theory as an input.
    If the overlap of the perturbed eigenstate with its corresponding unperturbed state are too small,
    this function raises an error, as perturbation theory breaks down.
    In this case, it also prints all states with a relevant admixture that should therefore be also included in the
    model space, to allow perturbation theory.
    """
    if not 0 < required_overlap <= 1:
        raise ValueError("Required overlap has to be a positive real number between zero and one.")

    model_inds = _get_model_inds(ket_tuples, system_pair)

    overlaps = (eigvec_perturb.multiply(eigvec_perturb.conj())).real  # elementwise multiplication
    for i, m_ind in enumerate(model_inds):
        overlaps[i, :] /= sparse.linalg.norm(overlaps[i, :])
        if overlaps[i, m_ind] >= required_overlap:
            continue
        logger.error(
            "The state %s has resonances with the following states, which might lead to unexpected/wrong results. "
            "Please consider adding them to your model space:",
            system_pair.basis.kets[m_ind],
        )
        print_above_admixture = (1 - required_overlap) * 0.05
        indices = sparse.find(overlaps[i, :] >= print_above_admixture)[1]
        for index in indices:
            if index == m_ind:
                continue
            admixture = 1 if np.isinf(overlaps[i, index]) else overlaps[i, index]
            logger.error("  - %s with admixture %.3f", system_pair.basis.kets[index], admixture)
