import logging
from collections.abc import Collection
from typing import TYPE_CHECKING, Any, Optional, Union, overload

import numpy as np
from scipy import sparse

from pairinteraction import (
    complex as pi_complex,
    real as pi_real,
)
from pairinteraction.units import QuantityArray, QuantityScalar, ureg

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from pint.facets.plain import PlainQuantity

    from pairinteraction._wrapped.ket.KetPair import KetPairLike

    KetAtom = Union[pi_real.KetAtom, pi_complex.KetAtom]
    KetPair = Union[pi_real.KetPair, pi_complex.KetPair]
    SystemPair = Union[pi_real.SystemPair, pi_complex.SystemPair]

logger = logging.getLogger(__name__)


@overload
def get_effective_hamiltonian(
    ket_tuple_list: Collection["KetPairLike"],
    magnetic_field: "PlainQuantity[NDArray[Any]]",
    electric_field: "PlainQuantity[NDArray[Any]]",
    distance_vector: "PlainQuantity[NDArray[Any]]",
    multipole_order: int = 3,
    perturbative_order: int = 2,
    with_diamagnetism: bool = False,
    required_overlap: float = 0.9,
    separate: bool = False,
) -> tuple["PlainQuantity[NDArray[Any]]", sparse.csr_matrix]: ...


@overload
def get_effective_hamiltonian(
    ket_tuple_list: Collection["KetPairLike"],
    magnetic_field: "PlainQuantity[NDArray[Any]]",
    electric_field: "PlainQuantity[NDArray[Any]]",
    distance_vector: "PlainQuantity[NDArray[Any]]",
    multipole_order: int = 3,
    perturbative_order: int = 2,
    with_diamagnetism: bool = False,
    required_overlap: float = 0.9,
    separate: bool = False,
    *,
    unit: str,
) -> tuple["NDArray[Any]", sparse.csr_matrix]: ...


def get_effective_hamiltonian(
    ket_tuple_list: Collection["KetPairLike"],
    magnetic_field: "PlainQuantity[NDArray[Any]]",
    electric_field: "PlainQuantity[NDArray[Any]]",
    distance_vector: "PlainQuantity[NDArray[Any]]",
    multipole_order: int = 3,
    perturbative_order: int = 2,
    with_diamagnetism: bool = False,
    required_overlap: float = 0.9,
    separate: bool = False,
    unit: Optional[str] = None,
) -> tuple[Union["NDArray[Any]", "PlainQuantity[NDArray[Any]]"], sparse.csr_matrix]:
    r"""Get the perturbative Hamiltonian at a desired order in Rayleigh-Schrödinger perturbation theory.

    This function takes a list of tuples of ket states, which forms the basis of the model space in which the effective
    Hamiltonian is calculated.
    The Hamiltonian of the pair system is assumed to be diagonal without atom-atom interactions, all atom-atom
    interaction is treated as a perturbative term.
    The function also checks for resonances between all states and states in the model space.

    Args:
        ket_tuple_list: List of all pair states that span up the model space.
            The effective Hamiltonian is calculated for these states.
        magnetic_field: magnetic field in the system.
        electric_field: electric field in the system.
        distance_vector: distance vector between the atoms.
        multipole_order: multipole-order of the interaction. Default is 3 (dipole-dipole).
        perturbative_order: order of perturbative calculation the system shall be used for. Default is 2.
        with_diamagnetism: True if diamagnetic term should be considered. Default is False.
        required_overlap: If set, the code checks for validity of a perturbative treatment.
            Error is thrown if the perturbed eigenstate has less overlap than this value with the
            unperturbed eigenstate.
        separate: True if only required order correction of the Hamiltonian is returned separately. False by default.
        unit: The unit to which to convert the result. Default None will return a pint quantity.

    Returns:
        - Effective Hamiltonian as a :math:`m \times m` matrix, where m is the length of `ket_tuple_list`.
        - Eigenvectors in perturbation theory due to interaction with states out of the model space, returned as
          a sparse matrix in compressed row format. Each row represents the corresponding eigenvector.

    Raises:
        ValueError: If a resonance between a state in the model space and a state not in the model space occurs.
        ValueError: If the perturbative order is not 1,2, or 3.

    """
    system_pair = _create_system(
        ket_tuple_list,
        magnetic_field,
        electric_field,
        distance_vector,
        multipole_order,
        perturbative_order,
        with_diamagnetism,
    )
    return get_effective_hamiltonian_from_system(
        ket_tuple_list,
        system_pair,
        order=perturbative_order,
        required_overlap=required_overlap,
        separate=separate,
        unit=unit,
    )


@overload
def get_effective_hamiltonian_from_system(
    ket_tuple_list: Collection["KetPairLike"],
    system_pair: "SystemPair",
    order: int = 2,
    required_overlap: float = 0.9,
    separate: bool = False,
) -> tuple["PlainQuantity[NDArray[Any]]", sparse.csr_matrix]: ...


@overload
def get_effective_hamiltonian_from_system(
    ket_tuple_list: Collection["KetPairLike"],
    system_pair: "SystemPair",
    order: int = 2,
    required_overlap: float = 0.9,
    separate: bool = False,
    *,
    unit: str,
) -> tuple["NDArray[Any]", sparse.csr_matrix]: ...


def get_effective_hamiltonian_from_system(
    ket_tuple_list: Collection["KetPairLike"],
    system_pair: "SystemPair",
    order: int = 2,
    required_overlap: float = 0.9,
    separate: bool = False,
    unit: Optional[str] = None,
) -> tuple[Union["NDArray[Any]", "PlainQuantity[NDArray[Any]]"], sparse.csr_matrix]:
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
        separate: True if only required order correction of the Hamiltonian is returned separately. False by default.
        unit: The unit to which to convert the result. Default None will return a pint quantity.

    Returns:
        - Effective Hamiltonian as a :math:`m \times m` matrix, where m is the length of `ket_tuple_list`.
        - Eigenvectors in perturbation theory due to interaction with states out of the model space, returned as
          a sparse matrix in compressed row format. Each row represents the corresponding eigenvector.

    Raises:
        ValueError: If a resonance between a state in the model space and a state not in the model space occurs.
        ValueError: If the order of perturbation theory is not 0,1,2 or 3.

    """
    if np.isinf(system_pair.get_distance().magnitude):
        raise ValueError(
            "Pair system is initialized without a distance. "
            "Please set a distance for calculating an effective Hamiltonian."
        )

    model_inds = _get_model_inds(ket_tuple_list, system_pair)
    h_au = system_pair.get_hamiltonian().to_base_units().magnitude  # Hamiltonian in atomic units
    h_eff_au, eigvec_perturb = _calculate_perturbative_hamiltonian(h_au, model_inds, order, separate)
    if not 0 <= required_overlap <= 1:
        raise ValueError("Required overlap has to be a positive real number between zero and one.")
    if required_overlap > 0:
        _check_for_resonances(model_inds, eigvec_perturb, system_pair, required_overlap)

    h_eff = QuantityArray.from_base_unit(h_eff_au, "ENERGY").to_pint_or_unit(unit)
    return h_eff, eigvec_perturb


@overload
def get_c3(
    ket_tuple_list: Collection["KetPairLike"],
    magnetic_field: "PlainQuantity[NDArray[Any]]",
    electric_field: "PlainQuantity[NDArray[Any]]",
    distance_vector: "PlainQuantity[NDArray[Any]]",
    multipole_order: int = 3,
    with_diamagnetism: bool = False,
) -> "PlainQuantity[float]": ...


@overload
def get_c3(
    ket_tuple_list: Collection["KetPairLike"],
    magnetic_field: "PlainQuantity[NDArray[Any]]",
    electric_field: "PlainQuantity[NDArray[Any]]",
    distance_vector: "PlainQuantity[NDArray[Any]]",
    multipole_order: int = 3,
    with_diamagnetism: bool = False,
    *,
    unit: str,
) -> float: ...


def get_c3(
    ket_tuple_list: Collection["KetPairLike"],
    magnetic_field: "PlainQuantity[NDArray[Any]]",
    electric_field: "PlainQuantity[NDArray[Any]]",
    distance_vector: "PlainQuantity[NDArray[Any]]",
    multipole_order: int = 3,
    with_diamagnetism: bool = False,
    unit: Optional[str] = None,
) -> Union[float, "PlainQuantity[float]"]:
    r"""Calculate the :math:`C_3` coefficient for a list of two 2-tuples of single atom ket states.

    This function calculates the :math:`C_3` coefficient in the desired unit. The input is a list of two 2-tuples of
    single atom ket states. We use the convention :math:`\Delta E = \frac{C_3}{r^3}`.

    Args:
        ket_tuple_list: The input as a list of tuples of two states [(a,b),(c,d)],
            the :math:`C_3` coefficient is calculated for (a,b)->(c,d).
            If there are not exactly two tuples in the list, a ValueError is raised.
        magnetic_field: magnetic field in the system.
        electric_field: electric field in the system.
        distance_vector: distance vector between the atoms.
        multipole_order: multipole-order of the interaction. Default is 3 (dipole-dipole).
        with_diamagnetism: True if diamagnetic term should be considered. Default is False.
        unit: The unit to which to convert the result. Default None will return a pint quantity.

    Returns:
        The :math:`C_3` coefficient with its unit.

    Raises:
        ValueError: If a list of not exactly two tuples of single atom states is given.

    """
    system_pair = _create_system(
        ket_tuple_list,
        magnetic_field,
        electric_field,
        distance_vector,
        multipole_order,
        perturbative_order=1,
        with_diamagnetism=with_diamagnetism,
    )
    return get_c3_from_system(ket_tuple_list, system_pair, unit)


@overload
def get_c3_from_system(
    ket_tuple_list: Collection["KetPairLike"],
    system_pair: "SystemPair",
) -> "PlainQuantity[float]": ...


@overload
def get_c3_from_system(ket_tuple_list: Collection["KetPairLike"], system_pair: "SystemPair", unit: str) -> float: ...


def get_c3_from_system(
    ket_tuple_list: Collection["KetPairLike"], system_pair: "SystemPair", unit: Optional[str] = None
) -> Union[float, "PlainQuantity[float]"]:
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
        c3 = get_c3_from_system(ket_tuple_list, system_pair, unit)
        system_pair.set_distance_vector(old_distance_vector)
        return c3

    h_eff, _ = get_effective_hamiltonian_from_system(ket_tuple_list, system_pair, order=1)
    c3 = h_eff[0, 1] * r**3
    return QuantityScalar.from_pint(c3, "C3").to_pint_or_unit(unit)


@overload
def get_c6(
    ket_tuple: "KetPairLike",
    magnetic_field: "PlainQuantity[NDArray[Any]]",
    electric_field: "PlainQuantity[NDArray[Any]]",
    distance_vector: "PlainQuantity[NDArray[Any]]",
    multipole_order: int = 3,
    with_diamagnetism: bool = False,
) -> "PlainQuantity[float]": ...


@overload
def get_c6(
    ket_tuple: "KetPairLike",
    magnetic_field: "PlainQuantity[NDArray[Any]]",
    electric_field: "PlainQuantity[NDArray[Any]]",
    distance_vector: "PlainQuantity[NDArray[Any]]",
    multipole_order: int = 3,
    with_diamagnetism: bool = False,
    *,
    unit: str,
) -> float: ...


def get_c6(
    ket_tuple: "KetPairLike",
    magnetic_field: "PlainQuantity[NDArray[Any]]",
    electric_field: "PlainQuantity[NDArray[Any]]",
    distance_vector: "PlainQuantity[NDArray[Any]]",
    multipole_order: int = 3,
    with_diamagnetism: bool = False,
    unit: Optional[str] = None,
) -> Union[float, "PlainQuantity[float]"]:
    r"""Calculate the :math:`C_6` coefficient for a given tuple of ket states.

    This function calculates the :math:`C_6` coefficient in the desired unit. The input is a 2-tuple of single atom ket
    states.

    Args:
        ket_tuple: The input is a tuple repeating the same single atom state in the format (a,a).
        If a tuple with not exactly two identical states is given, a ValueError is raised.
        magnetic_field: magnetic field in the system.
        electric_field: electric field in the system.
        distance_vector: distance vector between the atoms.
        multipole_order: multipole-order of the interaction. Default is 3 (dipole-dipole).
        with_diamagnetism: True if diamagnetic term should be considered. Default is False.
        unit: The unit to which to convert the result. Default None will return a pint quantity.

    Returns:
        The :math:`C_6` coefficient. If a unit is specified, the value in this unit is returned.

    Raises:
        ValueError: If a tuple with more than two single atom states is given.

    """
    system_pair = _create_system(
        [ket_tuple],
        magnetic_field,
        electric_field,
        distance_vector,
        multipole_order,
        perturbative_order=2,
        with_diamagnetism=with_diamagnetism,
    )
    return get_c6_from_system(ket_tuple, system_pair, unit=unit)


@overload
def get_c6_from_system(
    ket_tuple: "KetPairLike",
    system_pair: "SystemPair",
) -> "PlainQuantity[float]": ...


@overload
def get_c6_from_system(ket_tuple: "KetPairLike", system_pair: "SystemPair", unit: str) -> float: ...


def get_c6_from_system(
    ket_tuple: "KetPairLike", system_pair: "SystemPair", unit: Optional[str] = None
) -> Union[float, "PlainQuantity[float]"]:
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
        c6 = get_c6_from_system(ket_tuple, system_pair, unit)
        system_pair.set_distance_vector(old_distance_vector)
        return c6

    h_eff, _ = get_effective_hamiltonian_from_system([ket_tuple], system_pair, order=2, separate=True)
    c6 = h_eff[0, 0] * r**6
    return QuantityScalar.from_pint(c6, "C6").to_pint_or_unit(unit)


def _calculate_perturbative_hamiltonian(
    h: sparse.csr_matrix, model_inds: list[int], order: int = 2, separate: bool = False
) -> tuple["NDArray[Any]", sparse.csr_matrix]:
    r"""Calculate the perturbative Hamiltonian at a given order.

    This function takes a Hamiltonian as a sparse matrix which is diagonal in the unperturbed basis
    and list of indices spanning up the model space.
    It calculates both the effective Hamiltonian, spanned up by the states of the model space, as well as the
    perturbed eigenstates due to interactions with the exterior space in the desired order of perturbation theory.
    The output is either the full Hamiltonian up to the order of perturbation theory, or only the corrections at
    a given order.

    Args:
        h: Quadratic hermitian matrix. Perturbative terms are assumed to be only off-diagonal.
        model_inds: List of indices corresponding to the states that span up the model space.
        order: Order up to which the perturbation theory is expanded. Support up to third order.
            Default is second order.
        separate: True if a fixed order correction of the Hamiltonian is returned separately. False by default.

    Returns:
        Effective Hamiltonian as a :math:`m \times m` matrix, where m is the length of `ket_tuple_list`
        Eigenvectors in perturbation theory due to interaction with states out of the model
            space, returned as a sparse matrix in compressed row format. Each row represent the
            corresponding eigenvector

    """
    if order not in [0, 1, 2, 3]:
        raise ValueError("Perturbation theory is only implemented for orders [0, 1, 2, 3].")

    m_inds: NDArray[Any] = np.array(model_inds)
    o_inds: NDArray[Any] = np.setdiff1d(np.arange(h.shape[0]), m_inds)
    m = len(m_inds)
    n = len(o_inds)

    h0 = h.diagonal()
    h0_m = h0[m_inds]
    h_eff = [sparse.diags(h0_m, format="csr")]
    eigvec_perturb = sparse.hstack([sparse.eye(m, m, format="csr"), sparse.csr_matrix((m, n))])

    if order >= 1:
        v = h - sparse.diags(h0)
        v_mm = v[np.ix_(m_inds, m_inds)]
        h_eff.append(v_mm)

    if order >= 2:
        h0_e = h0[o_inds]
        v_me = v[np.ix_(m_inds, o_inds)]
        delta_energy_em = 1 / (h0_m[np.newaxis, :] - h0_e[:, np.newaxis])
        h_eff.append(v_me @ ((v_me.conj().T).multiply(delta_energy_em)))
        addition_mm = sparse.csr_matrix((m, m))
        addition_me = sparse.csr_matrix(((v_me.conj().T).multiply(delta_energy_em)).T)
        eigvec_perturb += sparse.hstack([addition_mm, addition_me])

    if order >= 3:
        diff = h0_m[np.newaxis, :] - h0_m[:, np.newaxis]
        diff = np.where(diff == 0, np.inf, diff)
        delta_energy_mm = 1 / diff
        v_ee = v[np.ix_(o_inds, o_inds)]
        if m > 1:
            logger.warning(
                "At third order, the eigenstates are currently only valid when only one state is in the model space. "
                "Take care with interpreation of the perturbed eigenvectors."
            )
        h_eff.append(
            v_me
            @ (
                (
                    v_ee @ ((v_me.conj().T).multiply(delta_energy_em))
                    - ((v_me.conj().T).multiply(delta_energy_em)) @ v_mm
                ).multiply(delta_energy_em)
            )
        )
        addition_mm_diag = -0.5 * sparse.diags(
            (v_me @ ((v_me.conj().T).multiply(np.square(delta_energy_em)))).diagonal(),
            format="csr",
        )
        addition_mm_offdiag = sparse.csr_matrix(
            ((v_me @ (v_me.conj().T).multiply(delta_energy_em)).multiply(delta_energy_mm)).T
        )
        addition_me = sparse.csr_matrix(
            ((v_ee @ ((v_me.conj().T).multiply(delta_energy_em))).multiply(delta_energy_em)).T
        )
        addition_me_2 = sparse.csr_matrix(
            ((v_me.conj().T @ ((v_mm.conj().T).multiply(delta_energy_mm))).multiply(delta_energy_em)).T
        )
        eigvec_perturb += sparse.hstack([addition_mm_diag + addition_mm_offdiag, addition_me + addition_me_2])

    # resort eigvec to original order
    all_inds = np.append(m_inds, o_inds)
    all_inds_positions = np.argsort(all_inds)
    eigvec_perturb = eigvec_perturb[:, all_inds_positions]
    if separate:
        return (0.5 * (h_eff[order] + h_eff[order].conj().T)).todense(), eigvec_perturb
    else:
        h_eff = np.cumsum(h_eff, axis=0)
        return (0.5 * (h_eff[-1] + h_eff[-1].conj().T)).todense(), eigvec_perturb


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
    eigvec_perturb: sparse.csr_matrix,
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
            logger.error(f"  - {system_pair.basis.kets[index]} with admixture {admixture:.3f}")
    if error_flag:
        raise ValueError(
            "Error. Perturbative Calculation not possible due to resonances. "
            "Add more states to the model space or adapt your required overlap."
        )


def _create_system(
    ket_tuple_list: Collection["KetPairLike"],
    magnetic_field: "PlainQuantity[NDArray[Any]]",
    electric_field: "PlainQuantity[NDArray[Any]]",
    distance_vector: "PlainQuantity[NDArray[Any]]",
    multipole_order: int = 3,
    perturbative_order: int = 2,
    with_diamagnetism: bool = False,
) -> "SystemPair":
    r"""Create a good estimate for a system in which to perform perturbation theory.

    This function takes a list of 2-tuples of ket states and creates a pair system holding a larger basis.
    The parameters of the basis are adjusted by the electric and magnetic field vectors of the system, as well
    as the distance and the multipole-order of the interaction. For higher-order perturbation theory, larger
    systems can be considered. Diamagnetism can be considered as well.

    Args:
        ket_tuple_list: List of all pair states that span up the model space. The system is created such that
        the effective Hamiltonian of the model system can be calculated accurately at a later stage.
        magnetic_field: magnetic field in the system.
        electric_field: electric field in the system.
        distance_vector: distance vector between the atoms.
        multipole_order: multipole-order of the interaction. Default is 3 (dipole-dipole).
        perturbative_order: order of perturbative calculation the system shall be used for. Default is 2.
        with_diamagnetism: True if diamagnetic term should be considered. Default is False.

    Returns:
        Pair system that can be used for perturbative calculations.

    Raises:
    ValueError: If the perturbative order is not 1,2 or 3.

    """
    if perturbative_order not in [1, 2, 3]:
        raise ValueError("System can only be created for perturbative order 1,2, and 3.")
    n_1 = []
    l_1 = []
    j_1 = []
    m_1 = []
    energ_au = []
    n_2 = []
    l_2 = []
    j_2 = []
    m_2 = []

    for ket1, ket2 in ket_tuple_list:
        n_1.append(ket1.n)
        l_1.append(ket1.l)
        j_1.append(ket1.j)
        m_1.append(ket1.m)
        n_2.append(ket2.n)
        l_2.append(ket2.l)
        j_2.append(ket2.j)
        m_2.append(ket2.m)
        energ_au.append(ket1.get_energy().magnitude + ket2.get_energy().magnitude)
    species = [ket_tuple_list[0][0].species, ket_tuple_list[0][1].species]
    spin = [ket_tuple_list[0][0].s, ket_tuple_list[0][1].s]
    n_max = [np.max(n_1), np.max(n_2)]
    l_max = [np.max(l_1), np.max(l_2)]
    j_max = [np.max(j_1), np.max(j_2)]
    m_max = [np.max(m_1), np.max(m_2)]
    n_min = [np.min(n_1), np.min(n_2)]
    l_min = [np.min(l_1), np.min(l_2)]
    j_min = [np.min(j_1), np.min(j_2)]
    m_min = [np.min(m_1), np.min(m_2)]
    ket_max = ket_tuple_list[np.argmax(energ_au)]
    ket_min = ket_tuple_list[np.argmin(energ_au)]

    if (
        magnetic_field[0].magnitude
        == magnetic_field[1].magnitude
        == electric_field[0].magnitude
        == electric_field[1].magnitude
        == 0
    ):
        pi = pi_real
        bases = [
            pi.BasisAtom(
                species[i],
                n=(n_min[i] - 7, n_max[i] + 7),
                l=(
                    np.maximum(0, l_min[i] - perturbative_order * (multipole_order - 2)),
                    l_max[i] + perturbative_order * (multipole_order - 2),
                ),
                j=(
                    np.maximum(0, j_min[i] - perturbative_order * (multipole_order - 2)),
                    j_max[i] + perturbative_order * (multipole_order - 2),
                ),
                m=(
                    m_min[i] - perturbative_order * (multipole_order - 2),
                    m_max[i] + perturbative_order * (multipole_order - 2),
                ),
            )
            for i in [0, 1]
        ]
    else:
        pi = pi_complex
        bases = [
            pi.BasisAtom(
                species[i],
                n=(n_min[i] - 7, n_max[i] + 7),
                l=(
                    np.maximum(0, l_min[i] - 1 - perturbative_order * (multipole_order - 2)),
                    l_max[i] + 1 + perturbative_order * (multipole_order - 2),
                ),
            )
            for i in [0, 1]
        ]
    distance_unit = distance_vector.units
    vector = distance_vector.magnitude
    distance = ureg.Quantity(np.linalg.norm(vector), distance_unit)
    population_admixture = 1e-4

    n_max_1 = ket_max[0].n
    n_max_2 = ket_max[1].n
    e_max = ket_max[0].get_energy() + ket_max[1].get_energy()
    dipole_max_1 = pi.KetAtom(
        species[0], n=n_max_1, l=n_max_1 - 1, j=n_max_1 - 1 + spin[0], m=n_max_1 - 1 + spin[0]
    ).get_matrix_element(
        pi.KetAtom(species[0], n=n_max_1 + 1, l=n_max_1, j=n_max_1 + spin[0], m=n_max_1 + spin[0]),
        operator="ELECTRIC_DIPOLE",
        q=1,
    )
    dipole_max_2 = pi.KetAtom(
        species[1], n=n_max_2, l=n_max_2 - 1, j=n_max_2 - 1 + spin[1], m=n_max_2 - 1 + spin[1]
    ).get_matrix_element(
        pi.KetAtom(species[1], n=n_max_2 + 1, l=n_max_2, j=n_max_2 + spin[1], m=n_max_2 + spin[1]),
        operator="ELECTRIC_DIPOLE",
        q=1,
    )
    shift_max = dipole_max_1 * dipole_max_2 * ureg.coulomb_constant / distance**3
    delta_energy_max = np.abs(shift_max) / np.sqrt(population_admixture) * perturbative_order
    n_min_1 = ket_min[0].n
    n_min_2 = ket_min[1].n
    e_min = ket_min[0].get_energy() + ket_min[1].get_energy()
    dipole_min_1 = pi.KetAtom(
        species[0], n=n_min_1, l=n_min_1 - 1, j=n_min_1 - 1 + spin[0], m=n_min_1 - 1 + spin[0]
    ).get_matrix_element(
        pi.KetAtom(species[0], n=n_min_1 + 1, l=n_min_1, j=n_min_1 + spin[0], m=n_min_1 + spin[0]),
        operator="ELECTRIC_DIPOLE",
        q=1,
    )
    dipole_min_2 = pi.KetAtom(
        species[1], n=n_min_2, l=n_min_2 - 1, j=n_min_2 - 1 + spin[1], m=n_min_2 - 1 + spin[1]
    ).get_matrix_element(
        pi.KetAtom(species[1], n=n_min_2 + 1, l=n_min_2, j=n_min_2 + spin[1], m=n_min_2 + spin[1]),
        operator="ELECTRIC_DIPOLE",
        q=1,
    )
    shift_min = dipole_min_1 * dipole_min_2 * ureg.coulomb_constant / distance**3
    delta_energy_min = np.abs(shift_min) / np.sqrt(population_admixture) * perturbative_order

    systems = [pi.SystemAtom(basis=basis) for basis in bases]
    for system in systems:
        system.enable_diamagnetism(with_diamagnetism)
        system.set_magnetic_field(magnetic_field)
        system.set_electric_field(electric_field)
    if not all(system.is_diagonal for system in systems):
        pi.diagonalize(systems, diagonalizer="eigen", sort_by_energy=False)
    basis_pair = pi.BasisPair(systems, energy=(e_min - delta_energy_min, e_max + delta_energy_max))
    system_pair = pi.SystemPair(basis_pair)
    system_pair.set_distance_vector(distance_vector)
    system_pair.set_order(multipole_order)
    return system_pair
