# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from collections.abc import Collection, Iterable
from typing import TYPE_CHECKING, Optional, Union, overload

import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix

from pairinteraction import (
    complex as pi_complex,
    real as pi_real,
)
from pairinteraction.units import AtomicUnits, QuantityArray, QuantityScalar, ureg

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix

    from pairinteraction._wrapped.ket.ket_atom import KetAtom  # noqa: F401  # required for sphinx for KetPairLike
    from pairinteraction._wrapped.ket.ket_pair import (  # noqa: F401  # required for sphinx for KetPairLike
        KetAtomTuple,
        KetPairComplex,
        KetPairLike,
        KetPairReal,
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
    return_only_specified_order: bool = False,
    unit: None = None,
) -> tuple["PintArray", "csr_matrix"]: ...


@overload
def get_effective_hamiltonian_from_system(
    ket_tuple_list: Collection["KetPairLike"],
    system_pair: "SystemPair",
    order: int = 2,
    required_overlap: float = 0.9,
    return_only_specified_order: bool = False,
    *,
    unit: str,
) -> tuple["NDArray", "csr_matrix"]: ...


def get_effective_hamiltonian_from_system(
    ket_tuple_list: Collection["KetPairLike"],
    system_pair: "SystemPair",
    order: int = 2,
    required_overlap: float = 0.9,
    return_only_specified_order: bool = False,
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
        return_only_specified_order: If True, the returned effective Hamiltonian will only contain the specified order.
            Default is False, which returns the sum of all orders up to the specified order.
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
    h_eff_au, eigvec_perturb = _calculate_perturbative_hamiltonian(h_au, model_inds, order, return_only_specified_order)
    if not 0 <= required_overlap <= 1:
        raise ValueError("Required overlap has to be a positive real number between zero and one.")
    if required_overlap > 0:
        _check_for_resonances(model_inds, eigvec_perturb, system_pair, required_overlap)

    h_eff = QuantityArray.convert_au_to_user(h_eff_au, "energy", unit)
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

    h_eff, _ = get_effective_hamiltonian_from_system(
        [ket_tuple], system_pair, order=2, return_only_specified_order=True
    )
    c6_pint = h_eff[0, 0] * r**6  # type: ignore [index] # PintArray does not know it can be indexed
    return QuantityScalar.from_pint(c6_pint, "c6").to_pint_or_unit(unit)


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


def create_system_for_perturbative(  # noqa: C901, PLR0912, PLR0915
    ket_tuple_list: Collection["KetAtomTuple"],
    electric_field: Optional["PintArray"] = None,
    magnetic_field: Optional["PintArray"] = None,
    distance_vector: Optional["PintArray"] = None,
    multipole_order: int = 3,
    with_diamagnetism: bool = False,
    perturbation_order: int = 2,
    number_of_considered_pair_kets: int = 2_000,
) -> "SystemPair":
    r"""Create a good estimate for a system in which to perform perturbative calculations.

    This function takes a list of 2-tuples of ket states and creates a pair system holding a larger basis.
    The parameters of the basis are adjusted by the electric and magnetic field vectors of the system, as well
    as the distance and the multipole-order of the interaction. For higher-order perturbation theory, larger
    systems can be considered. Diamagnetism can be considered as well.

    Args:
        ket_tuple_list: List of all pair states that span up the model space. The system is created such that
            the effective Hamiltonian of the model system can be calculated accurately at a later stage.
        electric_field: Electric field in the system.
        magnetic_field: Magnetic field in the system.
        distance_vector: Distance vector between the atoms.
        multipole_order: Multipole-order of the interaction. Default is 3 (dipole-dipole).
        with_diamagnetism: True if diamagnetic term should be considered. Default is False.
        perturbation_order: Order of perturbative calculation the system shall be used for. Default is 2.
        number_of_considered_pair_kets: Number of pair kets that are considered in the system. Default is 2000.

    Returns:
        Pair system that can be used for perturbative calculations.

    """
    electric_field = electric_field if electric_field is not None else ureg.Quantity([0, 0, 0], "V/cm")
    magnetic_field = magnetic_field if magnetic_field is not None else ureg.Quantity([0, 0, 0], "G")

    pi = pi_real if electric_field[1] == 0 and magnetic_field[1] == 0 else pi_complex  # type: ignore [index]
    are_fields_along_z = all(x == 0 for x in [*magnetic_field[:2], *electric_field[:2]])  # type: ignore [index]

    system_atoms: list[Union[pi_real.SystemAtom, pi_complex.SystemAtom]] = []

    delta_n = 7
    delta_l = perturbation_order * (multipole_order - 2)
    for i in range(2):
        kets = [ket_tuple[i] for ket_tuple in ket_tuple_list]
        nlfm = np.transpose([[ket.n, ket.l, ket.f, ket.m] for ket in kets])
        n_range = (int(np.min(nlfm[0])) - delta_n, int(np.max(nlfm[0])) + delta_n)
        l_range = (np.min(nlfm[1]) - delta_l, np.max(nlfm[1]) + delta_l)
        if any(ket.is_calculated_with_mqdt for ket in kets):
            # for mqdt we increase delta_l by 1 to take into account the variance ...
            l_range = (np.min(nlfm[1]) - delta_l - 1, np.max(nlfm[1]) + delta_l + 1)
        m_range = None
        if are_fields_along_z:
            m_range = (np.min(nlfm[3]) - delta_l, np.max(nlfm[3]) + delta_l)
        basis = pi.BasisAtom(kets[0].species, n=n_range, l=l_range, m=m_range)
        system = pi.SystemAtom(basis)
        system.set_diamagnetism_enabled(with_diamagnetism)
        system.set_magnetic_field(magnetic_field)
        system.set_electric_field(electric_field)
        system_atoms.append(system)

    pi.diagonalize(system_atoms)

    pair_energies_au = [
        sum(
            system.get_corresponding_energy(ket).to_base_units().magnitude
            for system, ket in zip(system_atoms, ket_tuple)
        )
        for ket_tuple in ket_tuple_list
    ]

    def get_basis_pair(delta_energy_au: float) -> Union[pi_real.BasisPair, pi_complex.BasisPair]:
        return pi.BasisPair(  # type: ignore [no-any-return]
            system_atoms,
            energy=(min(pair_energies_au) - delta_energy_au, max(pair_energies_au) + delta_energy_au),
            energy_unit=str(AtomicUnits["energy"]),
        )

    mhz_au = QuantityScalar.convert_user_to_au(1, "MHz", "energy")
    delta_energy_au = mhz_au
    min_delta, max_delta = None, None

    # make a bisect search to get a sensible basis size between:
    # number_of_considered_pair_kets and 1.2 * number_of_considered_pair_kets
    while delta_energy_au < 1:  # stop if delta_energy_au is 1 and the basis is still very small
        basis_pair = get_basis_pair(delta_energy_au)
        if basis_pair.number_of_kets < number_of_considered_pair_kets:
            min_delta = delta_energy_au
            if max_delta is None:
                delta_energy_au *= 2
            else:
                delta_energy_au = (delta_energy_au + max_delta) / 2
        elif basis_pair.number_of_kets > number_of_considered_pair_kets * 1.2:
            max_delta = delta_energy_au
            if min_delta is None:
                delta_energy_au /= 2
            else:
                delta_energy_au = (delta_energy_au + min_delta) / 2
        else:
            break
        if max_delta is not None and min_delta is not None and max_delta - min_delta < mhz_au:
            break

    logger.debug("The pair basis for the perturbative calculations consists of %d kets.", basis_pair.number_of_kets)

    system_pair = pi.SystemPair(basis_pair)
    if distance_vector is not None:
        system_pair.set_distance_vector(distance_vector)
    system_pair.set_interaction_order(multipole_order)
    return system_pair  # type: ignore [no-any-return]
