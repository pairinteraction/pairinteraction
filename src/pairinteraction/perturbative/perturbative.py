# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from collections.abc import Collection, Iterable
from typing import TYPE_CHECKING, Optional, Union, overload

import numpy as np
from scipy import sparse

from pairinteraction.units import QuantityArray, QuantityScalar

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix

    from pairinteraction import (
        complex as pi_complex,
        real as pi_real,
    )
    from pairinteraction._wrapped.ket.ket_atom import KetAtom  # noqa: F401  # required for sphinx for KetPairLike
    from pairinteraction._wrapped.ket.ket_pair import (
        KetPairComplex,  # noqa: F401  # required for sphinx for KetPairLike
        KetPairLike,
        KetPairReal,  # noqa: F401  # required for sphinx for KetPairLike
    )
    from pairinteraction.units import NDArray, PintArray, PintFloat

    SystemPair = Union[pi_real.SystemPair, pi_complex.SystemPair]

logger = logging.getLogger(__name__)


@overload
def get_effective_hamiltonian_from_system(
    ket_tuple_list: Collection["KetPairLike"],
    system_pair: "SystemPair",
    order: int = 2,
    required_overlap: float = 0.9,
    *,
    unit: None = None,
) -> tuple["PintArray", "csr_matrix"]: ...


@overload
def get_effective_hamiltonian_from_system(
    ket_tuple_list: Collection["KetPairLike"],
    system_pair: "SystemPair",
    order: int = 2,
    required_overlap: float = 0.9,
    *,
    unit: str,
) -> tuple["NDArray", "csr_matrix"]: ...


def get_effective_hamiltonian_from_system(
    ket_tuple_list: Collection["KetPairLike"],
    system_pair: "SystemPair",
    order: int = 2,
    required_overlap: float = 0.9,
    unit: Optional[str] = None,
) -> tuple[Union["NDArray", "PintArray"], "csr_matrix"]:
    r"""Get the perturbative Hamiltonian at a desired order in Rayleigh-Schrödinger perturbation theory.

    This function takes a list of tuples of ket states, which forms the basis of the model space in which the effective
    Hamiltonian is calculated. The whole Hamiltonian is taken from a pair system.
    The Hamiltonian of the pair system is assumed to be diagonal in the unperturbed Hamiltonian, all off-diagonal
    elements are assumed to belong to the perturbative term.
    The function also checks for resonances between all states and states in the model space.

    Args:
        ket_tuple_list: List of all pair states that span up the model space.
            The effective Hamiltonian is calculated for these states.
        system_pair: Two-Atom-System, diagonal in the basis of the unperturbed Hamiltonian.
        order: Order up to which the perturbation theory is expanded. Support up to third order.
            Default is second order.
        required_overlap: If set, the code checks for validity of a perturbative treatment.
            Error is thrown if the perturbed eigenstate has less overlap than this value with the
            unperturbed eigenstate.
        unit: The unit to which to convert the result. Default None will return a pint quantity.

    Returns:
        - Effective Hamiltonian as a :math:`m \times m` matrix, where m is the length of `ket_tuple_list`.
        - Eigenvectors in perturbation theory due to interaction with states out of the model space, returned as
          a sparse matrix in compressed row format. Each row represents the corresponding eigenvector.

    Raises:
        ValueError: If a resonance between a state in the model space and a state not in the model space occurs.

    """
    if np.isinf(system_pair.get_distance().magnitude):
        raise ValueError(
            "Pair system is initialized without a distance. "
            "Please set a distance for calculating an effective Hamiltonian."
        )

    model_inds = _get_model_inds(ket_tuple_list, system_pair)
    h_au = system_pair.get_hamiltonian().to_base_units().magnitude  # Hamiltonian in atomic units
    h_eff_au, eigvec_perturb = _calculate_perturbative_hamiltonian(h_au, model_inds, order)
    if not 0 <= required_overlap <= 1:
        raise ValueError("Required overlap has to be a positive real number between zero and one.")
    if required_overlap > 0:
        _check_for_resonances(model_inds, eigvec_perturb, system_pair, required_overlap)

    h_eff = QuantityArray.from_base_unit(h_eff_au, "energy").to_pint_or_unit(unit)
    return h_eff, eigvec_perturb


@overload
def get_c3_from_system(
    ket_tuple_list: Collection["KetPairLike"], system_pair: "SystemPair", *, unit: None = None
) -> "PintFloat": ...


@overload
def get_c3_from_system(ket_tuple_list: Collection["KetPairLike"], system_pair: "SystemPair", unit: str) -> float: ...


def get_c3_from_system(
    ket_tuple_list: Collection["KetPairLike"], system_pair: "SystemPair", unit: Optional[str] = None
) -> Union[float, "PintFloat"]:
    r"""Calculate the :math:`C_3` coefficient for a list of two 2-tuples of single atom ket states.

    This function calculates the :math:`C_3` coefficient in the desired unit. The input is a list of two 2-tuples of
    single atom ket states. We use the convention :math:`\Delta E = \frac{C_3}{r^3}`.

    Args:
        ket_tuple_list: The input as a list of tuples of two states [(a,b),(c,d)],
            the :math:`C_3` coefficient is calculated for (a,b)->(c,d).
            If there are not exactly two tuples in the list, a ValueError is raised.
        system_pair: The pair system that is used for the calculation.
        unit: The unit to which to convert the result. Default None will return a pint quantity.

    Returns:
        The :math:`C_3` coefficient with its unit.

    Raises:
        ValueError: If a list of not exactly two tuples of single atom states is given.

    """
    if len(ket_tuple_list) != 2:
        raise ValueError("C3 coefficient can be calculated only between two 2-atom states.")

    r = system_pair.get_distance()
    if np.isinf(r.magnitude):
        logger.warning(
            "Pair system is initialized without a distance. "
            "Calculating the C3 coefficient at a distance vector of [0, 0, 20] mum."
        )
        old_distance_vector = system_pair.get_distance_vector()
        system_pair.set_distance_vector([0, 0, 20], "micrometer")
        c3 = get_c3_from_system(ket_tuple_list, system_pair, unit=unit)
        system_pair.set_distance_vector(old_distance_vector)
        return c3

    h_eff, _ = get_effective_hamiltonian_from_system(ket_tuple_list, system_pair, order=1)
    c3_pint = h_eff[0, 1] * r**3  # type: ignore [index] # PintArray does not know it can be indexed
    return QuantityScalar.from_pint(c3_pint, "c3").to_pint_or_unit(unit)


@overload
def get_c6_from_system(ket_tuple: "KetPairLike", system_pair: "SystemPair", *, unit: None = None) -> "PintFloat": ...


@overload
def get_c6_from_system(ket_tuple: "KetPairLike", system_pair: "SystemPair", unit: str) -> float: ...


def get_c6_from_system(
    ket_tuple: "KetPairLike", system_pair: "SystemPair", unit: Optional[str] = None
) -> Union[float, "PintFloat"]:
    r"""Calculate the :math:`C_6` coefficient for a given tuple of ket states.

    This function calculates the :math:`C_6` coefficient in the desired unit. The input is a 2-tuple of single atom ket
    states.

    Args:
        ket_tuple: The input is a tuple repeating the same single atom state in the format (a,a).
        If a tuple with not exactly two identical states is given, a ValueError is raised.
        system_pair: The pair system that is used for the calculation.
        unit: The unit to which to convert the result. Default None will return a pint quantity.

    Returns:
        The :math:`C_6` coefficient. If a unit is specified, the value in this unit is returned.

    Raises:
        ValueError: If a tuple with more than two single atom states is given.

    """
    if isinstance(ket_tuple, Iterable):
        if len(ket_tuple) != 2:
            raise ValueError("C6 coefficient can be calculated only for a single 2-atom state.")
        if ket_tuple[0].species == ket_tuple[1].species and ket_tuple[0] != ket_tuple[1]:
            raise ValueError(
                "If you want to calculate 2nd order perturbations of two different states a and b, "
                "please use the get_effective_hamiltonian_from_system([(a,b), (b,a)], system_pair) function."
            )

    r = system_pair.get_distance()
    if np.isinf(r.magnitude):
        logger.warning(
            "Pair system is initialized without a distance. "
            "Calculating the C6 coefficient at a distance vector of [0, 0, 20] mum."
        )
        old_distance_vector = system_pair.get_distance_vector()
        system_pair.set_distance_vector([0, 0, 20], "micrometer")
        c6 = get_c6_from_system(ket_tuple, system_pair, unit=unit)
        system_pair.set_distance_vector(old_distance_vector)
        return c6

    h_eff, _ = get_effective_hamiltonian_from_system([ket_tuple], system_pair, order=2)
    h_0, _ = get_effective_hamiltonian_from_system([ket_tuple], system_pair, order=0)
    c6_pint = (h_eff[0, 0] - h_0[0, 0]) * r**6  # type: ignore [index] # PintArray does not know it can be indexed
    return QuantityScalar.from_pint(c6_pint, "c6").to_pint_or_unit(unit)


def _calculate_perturbative_hamiltonian(
    hamiltonian: "csr_matrix", model_inds: list[int], order: int = 2
) -> tuple["NDArray", "csr_matrix"]:
    r"""Calculate the perturbative Hamiltonian at a given order.

    This function takes a Hamiltonian as a sparse matrix which is diagonal in the unperturbed basis
    and list of indices spanning up the model space.
    It calculates both the effective Hamiltonian, spanned up by the states of the model space, as well as the
    perturbed eigenstates due to interactions with the exterior space in the desired order of perturbation theory.

    Args:
        hamiltonian: Quadratic hermitian matrix. Perturbative terms are assumed to be only off-diagonal.
        model_inds: List of indices corresponding to the states that span up the model space.
        order: Order up to which the perturbation theory is expanded. Support up to third order.
            Default is second order.

    Returns:
        Effective Hamiltonian as a :math:`m \times m` matrix, where m is the length of `ket_tuple_list`
        Eigenvectors in perturbation theory due to interaction with states out of the model
            space, returned as a sparse matrix in compressed row format. Each row represent the
            corresponding eigenvector

    """
    m_inds = np.array(model_inds)
    o_inds = np.setdiff1d(np.arange(hamiltonian.shape[0]), m_inds)
    h_eff, eigvec_perturb = _calculate_unsorted_perturbative_hamiltonian(hamiltonian, m_inds, o_inds, order)

    # resort eigvec to original order
    all_inds = np.append(m_inds, o_inds)
    all_inds_positions = np.argsort(all_inds)
    eigvec_perturb = eigvec_perturb[:, all_inds_positions]

    # include the hermitian conjugate part of the effective Hamiltonian
    h_eff = 0.5 * (h_eff + h_eff.conj().T)

    return h_eff, eigvec_perturb


def _calculate_unsorted_perturbative_hamiltonian(
    hamiltonian: "csr_matrix", m_inds: "NDArray", o_inds: "NDArray", order: int
) -> tuple["NDArray", "csr_matrix"]:
    # This function is outsourced from _calculate_perturbative_hamiltonian to allow for better type checking
    if order not in [0, 1, 2, 3]:
        raise ValueError("Perturbation theory is only implemented for orders [0, 1, 2, 3].")

    h0 = hamiltonian.diagonal()
    h0_m = h0[m_inds]

    h_eff = np.diag(h0_m)
    eigvec_perturb = sparse.hstack(
        [sparse.eye(len(m_inds), len(m_inds), format="csr").tocsr(), sparse.csr_matrix((len(m_inds), len(o_inds)))]
    ).tocsr()

    if order < 1:
        return h_eff, eigvec_perturb

    v_offdiag = hamiltonian - sparse.diags(h0)
    v_mm = v_offdiag[np.ix_(m_inds, m_inds)]
    h_eff += v_mm

    if order < 2:
        return h_eff, eigvec_perturb

    h0_e = h0[o_inds]
    v_me = v_offdiag[np.ix_(m_inds, o_inds)]
    delta_e_em = 1 / (h0_m[np.newaxis, :] - h0_e[:, np.newaxis])
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

    return h_eff, eigvec_perturb


def _get_model_inds(ket_tuple_list: Collection["KetPairLike"], system_pair: "SystemPair") -> list[int]:
    """Get the indices of all ket tuples in the basis of pair system.

    This function takes a list of 2-tuples of ket states, and a pair system holding the entire basis.
    It returns an array of indices of the pair system basis in the order of the tuple list.

    Args:
        ket_tuple_list: List of all pair states that span up the model space.
        system_pair: Two-Atom-System, diagonal in the basis of the unperturbed Hamiltonian.

    Returns:
        List of indices corresponding to the states that span up the model space.

    """
    model_inds = []
    for kets in ket_tuple_list:
        overlap = system_pair.basis.get_overlaps(kets)
        index = np.argmax(overlap)
        if overlap[index] == 0:
            raise ValueError(f"The pairstate {kets} is not part of the basis of the pair system.")
        if overlap[index] < 0.5:
            raise ValueError(f"The pairstate {kets} cannot be identified uniquely (max overlap: {overlap[index]}).")
        model_inds.append(int(index))
    return model_inds


def _check_for_resonances(
    model_inds: list[int],
    eigvec_perturb: "csr_matrix",
    system_pair: "SystemPair",
    required_overlap: float,
) -> None:
    r"""Check for resonance between the states in the model space and other states.

    This function takes the perturbed eigenvectors of the perturbation theory as an input.
    If the overlap of the perturbed eigenstate with its corresponding unperturbed state are too small,
    this function raises an error, as perturbation theory breaks down.
    In this case, it also prints all states with a relevant admixture that should therefore be also included in the
    model space, to allow perturbation theory.

    Args:
        model_inds: List of indices corresponding to the states that span up the model space.
        eigvec_perturb: Sparse representation of the perturbed eigenstates in the desired order of
            perturbation theory. Each row corresponds to the eigestate according to `state model indices.`
        system_pair: Two-Atom-System, diagonal in the basis of the unperturbed Hamiltonian.
        order: Order up to which the perturbation theory is expanded. Support up to third order.
            Default is second order.
        required_overlap: If set, the code checks for validity of a perturbative treatment.
            Error is thrown if the perturbed eigenstate has less overlap than this value with the unperturbed eigenstate

    Returns:
        Effective Hamiltonian as a :math:`m \times m` matrix, where m is the length of `ket_tuple_list`
        Eigenvectors in perturbation theory due to interaction with states out of the model space,
            returned as a sparse matrix in compressed row format. Each row represent the corresponding eigenvector

    Raises:
        ValueError: If a resonance between a state in the model space and a state not in the model space occurs.

    """
    overlaps = (eigvec_perturb.multiply(eigvec_perturb.conj())).real
    error_flag = False
    for i, j in zip(range(len(model_inds)), model_inds):
        vector_norm = sparse.linalg.norm(overlaps[i, :])
        overlap = overlaps[i, j] / vector_norm
        if overlap >= required_overlap:
            continue
        error_flag = True
        print_above_admixture = (1 - required_overlap) * 0.05
        indices = sparse.find(overlaps[i, :] >= print_above_admixture * vector_norm)[1]
        logger.error(
            "The state %s has resonances with the following states, please consider adding them to your model space:",
            system_pair.basis.kets[j],
        )
        for index in indices:
            if index == j:
                continue
            admixture = 1 if np.isinf(overlaps[i, index]) else overlaps[i, index] / vector_norm
            logger.error("  - %s with admixture %.3f", system_pair.basis.kets[index], admixture)
    if error_flag:
        raise ValueError(
            "Error. Perturbative Calculation not possible due to resonances. "
            "Add more states to the model space or adapt your required overlap."
        )
