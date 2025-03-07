import logging
from collections.abc import Iterable
from typing import Optional, Union

import numpy as np
import pint
from scipy import sparse

import pairinteraction.complex as pi_complex
import pairinteraction.real as pi_real
from pairinteraction.units import ureg

SystemPair = Union[pi_real.SystemPair, pi_complex.SystemPair]
KetAtom = Union[pi_real.KetAtom, pi_complex.KetAtom]


logger = logging.getLogger(__name__)


def get_effective_hamiltonian_from_system(
    ket_tuple_list: list[tuple[KetAtom, KetAtom]],
    system_pair: SystemPair,
    order: int = 2,
    required_overlap: Optional[float] = None,
    print_above_admixture: float = 1e-2,
    unit: Optional[str] = None,
) -> tuple[Union[np.ndarray, pint.Quantity], sparse.csr_matrix]:
    """Get the perturbative Hamiltonian at a desired order in Rayleigh-Schrödinger perturbation theory.

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
        unit: Unit in which the effective Hamiltonian is calculated. Default is GHz.
        required_overlap: If set, the code checks for validity of a perturbative treatment.
            Error is thrown if the perturbed eigenstate has less overlap than this value with the
            unperturbed eigenstate.
        print_above_admixture: If due to a resonance a state in the model space obtains an admixture larger than this
            value with a different basis state, the resonant states are logged.

    Returns:
        pint.Quantity: effective Hamiltonian as a mxm matrix with unit, where m is the length of `ket_tuple_list`.
        scipy.sparse.csr_matrix: eigenvectors in perturbation theory due to interaction with states out of the model
            space, returned as a sparse matrix in compressed row format. Each row represents the corresponding
            eigenvector.

    Raises:
        ValueError: If a resonance between a state in the model space and a state not in the model space occurs.

    """
    assert isinstance(ket_tuple_list, Iterable)
    model_space_indices = _get_model_space_indices(ket_tuple_list, system_pair)
    if unit:
        H = system_pair.get_hamiltonian(unit=unit)
    else:
        H = system_pair.get_hamiltonian().magnitude
    H_eff, eig_vec_perturb = _calculate_perturbative_hamiltonian(H, model_space_indices, order)
    if required_overlap:
        if required_overlap <= 0 or required_overlap > 1:
            raise ValueError("Required overlap has to be a positive real number between zero and one.")

        _check_for_resonances(
            model_space_indices, eig_vec_perturb, system_pair, required_overlap, print_above_admixture
        )
    if unit:
        return H_eff, eig_vec_perturb
    else:
        return ureg.Quantity(H_eff, units="bohr^2 * electron_mass/atomic_unit_of_time^2"), eig_vec_perturb


def get_c3_from_system(
    ket_tuple_list: list[tuple[KetAtom, KetAtom]], system_pair: SystemPair, unit: Optional[str] = None
) -> pint.Quantity:
    r"""Calculate the $C_3$ coefficient for a list of two 2-tuples of single atom ket states.

    This function calculates the $C_3$ coefficient in the desired unit. The input is a list of two 2-tuples of
    single atom ket states. We use the convention $\Delta E = \frac{C_3}{r^3}$.

    Args:
        ket_tuple_list: The input as a list of tuples of two states [(a,b),(c,d)], the $C_3$ coefficient is caluculated
        for (a,b)->(c,d).
        If there are not exactly two tuples in the list, a ValueError is raised.
        system_pair: The pair system that is used for the calculation.
        unit: Unit of of the $C_3$ coeffcient.

    Returns:
        ureg.Quantity: The $C_3$ coefficient with its unit. If a unit is specified, the value in this unit is returned.

    Raises:
        ValueError: If a list of not exactly two tuples of single atom states is given.

    """
    assert isinstance(ket_tuple_list, Iterable)
    if len(ket_tuple_list) != 2:
        raise ValueError("C3 coefficient can be calculated only between two 2-atom states.")
    # ket_tuple_list = list(itertools.product(ket_tuple_list, ket_tuple_list))[1:3]
    vector = np.array([comp.magnitude for comp in system_pair.distance_vector])
    r = ureg.Quantity(np.linalg.norm(vector), units="bohr")
    H_eff, _ = get_effective_hamiltonian_from_system(ket_tuple_list, system_pair, order=1)
    C3 = H_eff[0, 1] * r**3
    if unit:
        C3 = C3.to(unit)
    return C3


def get_c6_from_system(
    ket_tuple: tuple[KetAtom, KetAtom], system_pair: SystemPair, unit: Optional[str] = None
) -> pint.Quantity:
    r"""Calculate the $C_6$ coefficient for a given tuple of ket states.

    This function calculates the $C_6$ coefficient in the desired unit. The input is a 2-tuple of single atom ket
    states.

    Args:
        ket_tuple: The input is a tuple repeating the same single atom state in the format (a,a).
        If a tuple with not exactly two identical states is given, a ValueError is raised.
        system_pair: The pair system that is used for the calculation.
        unit: Desired unit of of the $C_6$ coeffcient.

    Returns:
        ureg.Quantity: The $C_6$ coefficient with its unit. If a unit is specified, the value in this unit is returned.

    Raises:
        ValueError: If a tuple with more than two single atom states is given.

    """
    assert isinstance(ket_tuple, Iterable)
    if len(ket_tuple) != 2:
        # @Johannes: Add a fast comparison for ket states that gives true for the same state.
        raise ValueError("C6 coefficient can be calculated only for a single 2-atom state.")
    vector = np.array([comp.magnitude for comp in system_pair.distance_vector])
    r = ureg.Quantity(np.linalg.norm(vector), units="bohr")
    H_eff, _ = get_effective_hamiltonian_from_system([ket_tuple], system_pair, order=2)
    H_0, _ = get_effective_hamiltonian_from_system([ket_tuple], system_pair, order=0)
    C6 = (H_eff - H_0)[0, 0] * r**6
    if unit:
        C6 = C6.to(unit)
    return C6


def _calculate_perturbative_hamiltonian(
    H: sparse.csr_matrix, model_space_indices: np.ndarray, order: int = 2
) -> tuple[np.ndarray, sparse.csr_matrix]:
    """Calculate the perturbative Hamiltonian at a given order.

    This function takes a Hamiltonian as a sparse matrix which is diagonal in the unperturbed basis
    and list of indices spanning up the model space.
    It calculates both the effective Hamiltonian, spanned up by the states of the model space, as well as the perturbed
    eigenstates due to interactions with the exterior space in the desired order of perturbation theory.

    Args:
        H: Quadratic hermitian matrix. Perturbative terms are assumed to be only off-diagonal.
        model_space_indices: List of indices corresponding to the states that span up the model space.
        order: Order up to which the perturbation theory is expanded. Support up to third order.
        Default is second order.

    Returns:
        np.ndarray: effective Hamiltonian as a mxm matrix, where m is the length of `ket_tuple_list`
        scipy.sparse.csr_matrix: eigenvectors in perturbation theory due to interaction with states out of the model
        space, returned as a sparse matrix in compressed row format. Each row represent the corresponding eigenvector.

    """
    if order < 0 or order > 3 or int(order) != order:
        raise ValueError("Perturbation theory is only implemented up to the third order.")
    other_indices = np.setdiff1d(np.arange(H.shape[0]), model_space_indices)
    m = len(model_space_indices)
    n = len(other_indices)
    all_indices = np.append(model_space_indices, other_indices)
    all_indices_positions = np.argsort(all_indices)

    H0 = H.diagonal()
    H0_m = H0[model_space_indices]
    H_eff = sparse.diags(H0_m, format="csr")
    eig_vec_perturb = sparse.hstack([sparse.eye(m, m, format="csr"), sparse.csr_matrix((m, n))])[
        :, all_indices_positions
    ]

    if order >= 1:
        V = H - sparse.diags(H0)
        V_mm = V[np.ix_(model_space_indices, model_space_indices)]
        H_eff += V_mm

    if order >= 2:
        H0_e = H0[other_indices]
        V_me = V[np.ix_(model_space_indices, other_indices)]
        deltaE_em = 1 / (H0_m[np.newaxis, :] - H0_e[:, np.newaxis])
        H_eff += V_me @ ((V_me.conj().T).multiply(deltaE_em))
        addition_mm = sparse.csr_matrix((m, m))
        addition_me = sparse.csr_matrix(((V_me.conj().T).multiply(deltaE_em)).T)
        eig_vec_perturb += sparse.hstack([addition_mm, addition_me])[:, all_indices_positions]

    if order >= 3:
        diff = H0_m[np.newaxis, :] - H0_m[:, np.newaxis]
        diff = np.where(diff == 0, np.inf, diff)
        deltaE_mm = 1 / diff
        V_ee = V[np.ix_(other_indices, other_indices)]
        if m > 1:
            logger.warning(
                "At third order, the eigenstates are currently only valid when only one state is in the model space."
                "Take care with interpreation of the perturbed eigenvectors."
            )
        H_eff += V_me @ (
            (V_ee @ ((V_me.conj().T).multiply(deltaE_em)) - ((V_me.conj().T).multiply(deltaE_em)) @ V_mm).multiply(
                deltaE_em
            )
        )
        addition_mm_diag = -0.5 * sparse.diags(
            (V_me @ ((V_me.conj().T).multiply(np.square(deltaE_em)))).diagonal(),
            format="csr",
        )
        addition_mm_offdiag = sparse.csr_matrix(((V_me @ (V_me.conj().T).multiply(deltaE_em)).multiply(deltaE_mm)).T)
        addition_me = sparse.csr_matrix(((V_ee @ ((V_me.conj().T).multiply(deltaE_em))).multiply(deltaE_em)).T)
        addition_me_2 = sparse.csr_matrix(
            ((V_me.conj().T @ ((V_mm.conj().T).multiply(deltaE_mm))).multiply(deltaE_em)).T
        )
        eig_vec_perturb += sparse.hstack([addition_mm_diag + addition_mm_offdiag, addition_me + addition_me_2])[
            :, all_indices_positions
        ]

    return (0.5 * (H_eff + H_eff.conj().T)).todense(), eig_vec_perturb


def _get_model_space_indices(ket_tuple_list: list[tuple[KetAtom, KetAtom]], system_pair: SystemPair) -> np.ndarray:
    """Get the indices of all ket tuples in the basis of pair system.

    This function takes a list of 2-tuples of ket states, and a pair system holding the entire basis.
    It returns an array of indices of the pair system basis in the order of the tuple list.

    Args:
        ket_tuple_list: List of all pair states that span up the model space.
        system_pair: Two-Atom-System, diagonal in the basis of the unperturbed Hamiltonian.

    Returns:
        np.ndarray: List of indices corresponding to the states that span up the model space.

    """
    model_space_indices = []
    for kets in ket_tuple_list:
        overlap = system_pair.basis.get_overlaps(kets)
        index = np.argmax(overlap)
        if overlap[index] < 0.5:
            raise ValueError("The pairstate |" + str(kets[0]) + ">|" + str(kets[1]) + "> cannot be identified uniquely")
        model_space_indices.append(index)
    return np.array(model_space_indices)


def _check_for_resonances(
    model_space_indices: np.ndarray,
    eig_vec_perturb: sparse.csr_matrix,
    system_pair: SystemPair,
    required_overlap: float,
    print_above_admixture: float,
) -> None:
    """Check for resonance between the states in the model space and other states.

    This function takes the perturbed eigenvectors of the perturbation theory as an input.
    If the overlap of the perturbed eigenstate with its corresponding unperturbed state are too small,
    this function raises an error, as perturbation theory breaks down.
    In this case, it also prints all states with a relevant admixture that should therefore be also included in the
    model space, to allow perturbation theory.

    Args:
        model_space_indices: List of indices corresponding to the states that span up the model space.
        eig_vec_perturb: Sparse representation of the perturbed eigenstates in the desired order of
        perturbation theory. Each row corresponds to the eigestate according to `state model indices.`
        system_pair: Two-Atom-System, diagonal in the basis of the unperturbed Hamiltonian.
        order: Order up to which the perturbation theory is expanded. Support up to third order.
        Default is second order.
        required_overlap: If set, the code checks for validity of a perturbative treatment.
        Error is thrown if the perturbed eigenstate has less overlap than this value with the unperturbed eigenstate.
        print_above_admixture: If due to a resonance a state in the model space obtains an admixture larger than
        this value with a different basis state, the resonant states are logged.

    Returns:
        np.ndarray: effective Hamiltonian as a mxm matrix, where m is the length of `ket_tuple_list`
        scipy.sparse.csr_matrix: eigenvectors in perturbation theory due to interaction with states out of the model
        space, returned as a sparse matrix in compressed row format. Each row represent the corresponding eigenvector.

    Raises:
        ValueError: If a resonance between a state in the model space and a state not in the model space occurs.

    """
    overlaps = (eig_vec_perturb.multiply(eig_vec_perturb.conj())).real
    for i, j in zip(np.arange(len(model_space_indices)), model_space_indices):
        vector_norm = sparse.linalg.norm(overlaps[i])
        overlap = overlaps[i, j] / vector_norm
        if overlap < required_overlap:
            error_flag = True
            logger.error(
                "Error. Perturbative Calculation not possible due to resonance of state |"
                + str(system_pair.basis.kets[j])
                + ">."
            )
            indices = sparse.find(overlaps[i] >= print_above_admixture * vector_norm)[1]
            for index in indices:
                if index != j:
                    logger.error(
                        "The error occurs due to a resonance between the states |"
                        + str(system_pair.basis.kets[j])
                        + "> and |"
                        + str(system_pair.basis.kets[index])
                        + ">. Please include |"
                        + str(system_pair.basis.kets[index])
                        + "> also in your model space."
                    )
    if error_flag:
        raise ValueError(
            "Perturbation Theory impossible due to resonances with states not encountered in the model space."
            " See above for further details."
        )
