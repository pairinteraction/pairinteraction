# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import copy
import logging
from collections.abc import Sequence
from functools import cached_property, lru_cache
from typing import TYPE_CHECKING, Any, Literal, Optional, Union, overload

import numpy as np
from scipy import sparse

import pairinteraction.real as pi_real
from pairinteraction.perturbative.perturbation_theory import calculate_perturbative_hamiltonian
from pairinteraction.units import QuantityArray, QuantityScalar

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix
    from typing_extensions import Self

    from pairinteraction._wrapped.basis.basis_atom import BasisAtom
    from pairinteraction._wrapped.basis.basis_pair import BasisPair
    from pairinteraction._wrapped.ket.ket_atom import (
        KetAtom,  # noqa: F401  # needed for sphinx to recognize KetAtomTuple
    )
    from pairinteraction._wrapped.ket.ket_pair import KetAtomTuple
    from pairinteraction._wrapped.system.system_atom import SystemAtom
    from pairinteraction._wrapped.system.system_pair import SystemPair
    from pairinteraction.units import ArrayLike, NDArray, PintArray, PintFloat


logger = logging.getLogger(__name__)


class EffectiveSystemPair:
    def __init__(self, ket_tuples: Sequence["KetAtomTuple"]) -> None:
        if not all(len(ket_tuple) == 2 for ket_tuple in ket_tuples):
            raise ValueError("All ket tuples must contain exactly two kets")
        for i in range(2):
            if not all(ket_tuple[i].species == ket_tuples[0][i].species for ket_tuple in ket_tuples):
                raise ValueError(f"All kets for atom={i} must have the same species")

        # Perturbation attributes
        self._ket_tuples = [tuple(kets) for kets in ket_tuples]
        self._perturbation_order = 2

        # BasisAtom and SystemAtom attributes
        self._electric_field: Optional[PintArray] = None
        self._magnetic_field: Optional[PintArray] = None
        self._diamagnetism_enabled: Optional[bool] = None

        # BasisPair and SystemPair attributes
        self._minimum_number_of_ket_pairs: Optional[int] = None
        self._maximum_number_of_ket_pairs: Optional[int] = None
        self._interaction_order: Optional[int] = None
        self._distance_vector: Optional[PintArray] = None

        # misc
        self._eff_h_dict_au: Optional[dict[int, NDArray]] = None
        self._eff_vecs: Optional[csr_matrix] = None

    def copy(self: "Self") -> "Self":
        """Create a copy of the EffectiveSystemPair object (before it has been created)."""
        if hasattr(self, "_basis_atoms"):
            raise RuntimeError(
                "Cannot copy the EffectiveSystemPair object after it has been created. "
                "Please create a new object instead."
            )
        return copy.copy(self)

    def _ensure_not_created(
        self: "Self", what: Literal["basis_atoms", "system_atoms", "basis_pair", "system_pair"] = "basis_atoms"
    ) -> None:
        """Check if the effective Hamiltonian has already been created."""
        if hasattr(self, "_" + what):
            raise RuntimeError(
                f"Cannot change parameters for {what} after it has already been created. "
                f"Please set all parameters before {what} before accessing it (or creating the effective Hamiltonian)."
            )

    # # # Perturbation methods and attributes # # #
    @property
    def ket_tuples(self) -> list["KetAtomTuple"]:
        """The tuples of kets, which form the model space for the effective Hamiltonian."""
        return self._ket_tuples  # type: ignore [return-value]

    @property
    def perturbation_order(self) -> int:
        """The perturbation order for the effective Hamiltonian."""
        return self._perturbation_order

    def set_perturbation_order(self: "Self", order: int) -> "Self":
        """Set the perturbation order for the effective Hamiltonian."""
        self._ensure_not_created()
        self._perturbation_order = order
        return self

    # # # BasisAtom methods and attributes # # #
    @property
    def basis_atoms(self) -> tuple["BasisAtom[Any]", "BasisAtom[Any]"]:
        """The basis objects for the single-atom systems."""
        if not hasattr(self, "_basis_atoms"):
            self._create_basis_atoms()
        return self._basis_atoms  # type: ignore [return-value]

    @basis_atoms.setter
    def basis_atoms(self, basis_atoms: tuple["BasisAtom[Any]", "BasisAtom[Any]"]) -> None:
        self._ensure_not_created("basis_atoms")
        self._basis_atoms = tuple(basis_atoms)

    def _create_basis_atoms(self) -> None:
        delta_n = 7
        delta_l = self.perturbation_order * (self.interaction_order - 2)

        basis_atoms: list[BasisAtom[Any]] = []
        for i in range(2):
            kets = [ket_tuple[i] for ket_tuple in self.ket_tuples]
            nlfm = np.transpose([[ket.n, ket.l, ket.f, ket.m] for ket in kets])
            n_range = (int(np.min(nlfm[0])) - delta_n, int(np.max(nlfm[0])) + delta_n)
            l_range = (np.min(nlfm[1]) - delta_l, np.max(nlfm[1]) + delta_l)
            if any(ket.is_calculated_with_mqdt for ket in kets):
                # for mqdt we increase delta_l by 1 to take into account the variance ...
                l_range = (np.min(nlfm[1]) - delta_l - 1, np.max(nlfm[1]) + delta_l + 1)
            m_range = None
            if self._are_fields_along_z:
                m_range = (np.min(nlfm[3]) - delta_l, np.max(nlfm[3]) + delta_l)
            basis = get_basis_atom_with_cache(kets[0].species, n_range, l_range, m_range)
            basis_atoms.append(basis)

        self._basis_atoms = tuple(basis_atoms)

    # # # SystemAtom methods and attributes # # #
    @property
    def system_atoms(self) -> tuple["SystemAtom[Any]", "SystemAtom[Any]"]:
        """The system objects for the single-atom systems."""
        if not hasattr(self, "_system_atoms"):
            self._create_system_atoms()
        return self._system_atoms  # type: ignore [return-value]

    @system_atoms.setter
    def system_atoms(self, system_atoms: tuple["SystemAtom[Any]", "SystemAtom[Any]"]) -> None:
        self._ensure_not_created()
        self._system_atoms = tuple(system_atoms)
        self.basis_atoms = tuple(system.basis for system in system_atoms)

    @property
    def electric_field(self) -> "PintArray":
        """The electric field for the single-atom systems."""
        if self._electric_field is None:
            self.set_electric_field([0, 0, 0], "V/cm")
            assert self._electric_field is not None
        return self._electric_field

    def set_electric_field(
        self: "Self",
        electric_field: Union["PintArray", "ArrayLike"],
        unit: Optional[str] = None,
    ) -> "Self":
        """Set the electric field for the single-atom systems.

        Args:
            electric_field: The electric field to set for the systems.
            unit: The unit of the electric field, e.g. "V/cm".
                Default None expects a `pint.Quantity`.

        """
        self._ensure_not_created()
        self._electric_field = QuantityArray.convert_user_to_pint(electric_field, unit, "electric_field")
        return self

    @property
    def magnetic_field(self) -> "PintArray":
        """The magnetic field for the single-atom systems."""
        if self._magnetic_field is None:
            self.set_magnetic_field([0, 0, 0], "gauss")
            assert self._magnetic_field is not None
        return self._magnetic_field

    def set_magnetic_field(
        self: "Self",
        magnetic_field: Union["PintArray", "ArrayLike"],
        unit: Optional[str] = None,
    ) -> "Self":
        """Set the magnetic field for the single-atom systems.

        Args:
            magnetic_field: The magnetic field to set for the systems.
            unit: The unit of the magnetic field, e.g. "gauss".
                Default None expects a `pint.Quantity`.

        """
        self._ensure_not_created()
        self._magnetic_field = QuantityArray.convert_user_to_pint(magnetic_field, unit, "magnetic_field")
        return self

    @property
    def _are_fields_along_z(self) -> bool:
        return all(x == 0 for x in [*self.magnetic_field[:2], *self.electric_field[:2]])  # type: ignore [index]

    @property
    def diamagnetism_enabled(self) -> bool:
        """Whether diamagnetism is enabled for the single-atom systems."""
        if self._diamagnetism_enabled is None:
            self.set_diamagnetism_enabled(False)
            assert self._diamagnetism_enabled is not None
        return self._diamagnetism_enabled

    def set_diamagnetism_enabled(self: "Self", enable: bool = True) -> "Self":
        """Enable or disable diamagnetism for the system.

        Args:
            enable: Whether to enable or disable diamagnetism.

        """
        self._ensure_not_created("system_atoms")
        self._diamagnetism_enabled = enable
        return self

    def _create_system_atoms(self) -> None:
        system_atoms: list[SystemAtom[Any]] = []
        for basis_atom in self.basis_atoms:
            system = pi_real.SystemAtom(basis_atom)  # type: ignore [arg-type]
            system.set_diamagnetism_enabled(self.diamagnetism_enabled)
            system.set_electric_field(self.electric_field)
            system.set_magnetic_field(self.magnetic_field)
            system_atoms.append(system)
        pi_real.diagonalize(system_atoms)

        self._system_atoms = tuple(system_atoms)

    @property
    def pair_energies(self) -> list["PintFloat"]:
        return [
            sum(system.get_corresponding_energy(ket) for system, ket in zip(self.system_atoms, ket_tuple, strict=True))  # type: ignore [misc]
            for ket_tuple in self.ket_tuples
        ]

    # # # BasisPair methods and attributes # # #
    @property
    def basis_pair(self) -> "BasisPair[Any, Any]":
        """The basis pair object for the pair system."""
        if not hasattr(self, "_basis_pair"):
            self._create_basis_pair()
        return self._basis_pair

    @basis_pair.setter
    def basis_pair(self, basis_pair: "BasisPair[Any, Any]") -> None:
        self._ensure_not_created()
        self._basis_pair = basis_pair
        self.system_atoms = basis_pair.system_atoms

    def set_minimum_number_of_ket_pairs(self: "Self", number_of_kets: int) -> "Self":
        """Set the minimum number of ket pairs in the basis pair.

        Args:
            number_of_kets: The minimum number of ket pairs to set in the basis pair, by default we use 2000.

        """
        self._ensure_not_created("basis_pair")
        if self._maximum_number_of_ket_pairs is not None and number_of_kets > self._maximum_number_of_ket_pairs:
            raise ValueError("The minimum number of ket pairs cannot be larger than the maximum number of ket pairs.")
        self._minimum_number_of_ket_pairs = number_of_kets
        return self

    def set_maximum_number_of_ket_pairs(self: "Self", number_of_kets: int) -> "Self":
        """Set the maximum number of ket pairs in the basis pair.

        Args:
            number_of_kets: The maximum number of ket pairs to set in the basis pair.

        """
        self._ensure_not_created("basis_pair")
        if self._minimum_number_of_ket_pairs is not None and number_of_kets < self._minimum_number_of_ket_pairs:
            raise ValueError("The maximum number of ket pairs cannot be smaller than the minimum number of ket pairs.")
        self._maximum_number_of_ket_pairs = number_of_kets
        return self

    def _create_basis_pair(self) -> None:
        max_number_of_kets: Optional[float] = self._maximum_number_of_ket_pairs
        if max_number_of_kets is None:
            max_number_of_kets = np.inf
        min_number_of_kets: Optional[float] = self._minimum_number_of_ket_pairs
        if min_number_of_kets is None:
            min_number_of_kets = min(2_000, max_number_of_kets)

        pair_energies_au = [energy.to_base_units().magnitude for energy in self.pair_energies]
        min_energy_au = min(pair_energies_au)
        max_energy_au = max(pair_energies_au)

        mhz_au = QuantityScalar.convert_user_to_au(1, "MHz", "energy")
        delta_energy_au = 100 * mhz_au
        min_delta, max_delta = 0.1 * mhz_au, 1.0

        # make a bisect search to get a sensible basis size between:
        # min_number_of_kets and max_number_of_kets
        while True:
            basis_pair = pi_real.BasisPair(
                self.system_atoms,
                energy=(min_energy_au - delta_energy_au, max_energy_au + delta_energy_au),
                energy_unit="hartree",
            )

            if max_delta - min_delta < mhz_au:
                break  # stop condition if delta_energy_au does not change anymore

            if basis_pair.number_of_kets < min_number_of_kets:
                min_delta = delta_energy_au
                delta_energy_au = min(2 * delta_energy_au, (delta_energy_au + max_delta) / 2)
            elif basis_pair.number_of_kets > max_number_of_kets:
                max_delta = delta_energy_au
                min_delta = max(0.5 * delta_energy_au, (delta_energy_au + min_delta) / 2)
            else:
                break

        self._basis_pair = basis_pair
        logger.debug("The pair basis for the perturbative calculations consists of %d kets.", basis_pair.number_of_kets)

    # # # SystemPair methods and attributes # # #
    @property
    def system_pair(self) -> "SystemPair[Any]":
        """The system pair object for the pair system."""
        if not hasattr(self, "_system_pair"):
            self._create_system_pair()
        return self._system_pair

    @system_pair.setter
    def system_pair(self, system_pair: "SystemPair[Any]") -> None:
        self._ensure_not_created()
        self._system_pair = system_pair
        self.basis_pair = system_pair.basis

    @property
    def interaction_order(self) -> int:
        """The interaction order for the pair system."""
        if self._interaction_order is None:
            self.set_interaction_order(3)
        return self._interaction_order  # type: ignore [return-value]

    def set_interaction_order(self: "Self", order: int) -> "Self":
        """Set the interaction order of the pair system.

        Args:
            order: The interaction order to set for the pair system.
                The order must be 3, 4, or 5.

        """
        self._ensure_not_created()
        self._interaction_order = order
        return self

    @property
    def distance_vector(self) -> "PintArray":
        """The distance vector between the atoms in the pair system."""
        if self._distance_vector is None:
            self.set_distance_vector([0, 0, np.inf], "micrometer")
        return self._distance_vector  # type: ignore [return-value]

    def set_distance(
        self: "Self",
        distance: Union[float, "PintFloat"],
        angle_degree: float = 0,
        unit: Optional[str] = None,
    ) -> "Self":
        """Set the distance between the atoms using the specified distance and angle.

        Args:
            distance: The distance to set between the atoms in the given unit.
            angle_degree: The angle between the distance vector and the z-axis in degrees.
                90 degrees corresponds to the x-axis.
                Defaults to 0, which corresponds to the z-axis.
            unit: The unit of the distance, e.g. "micrometer".
                Default None expects a `pint.Quantity`.

        """
        distance_vector = [np.sin(np.deg2rad(angle_degree)) * distance, 0, np.cos(np.deg2rad(angle_degree)) * distance]
        return self.set_distance_vector(distance_vector, unit)

    def set_distance_vector(
        self: "Self",
        distance: Union["ArrayLike", "PintArray"],
        unit: Optional[str] = None,
    ) -> "Self":
        """Set the distance vector between the atoms.

        Args:
            distance: The distance vector to set between the atoms in the given unit.
            unit: The unit of the distance, e.g. "micrometer".
                Default None expects a `pint.Quantity`.

        """
        self._ensure_not_created("system_pair")
        self._distance_vector = QuantityArray.convert_user_to_pint(distance, unit, "distance")
        return self

    def _create_system_pair(self) -> None:
        system_pair = pi_real.SystemPair(self.basis_pair)  # type: ignore [arg-type]
        system_pair.set_distance_vector(self.distance_vector)
        system_pair.set_interaction_order(self.interaction_order)
        self._system_pair = system_pair

    # # # Effective Hamiltonian methods and attributes # # #
    @overload
    def get_effective_hamiltonian(self, return_order: Optional[int] = None, unit: None = None) -> "PintArray": ...

    @overload
    def get_effective_hamiltonian(self, return_order: Optional[int] = None, *, unit: str) -> "NDArray": ...

    def get_effective_hamiltonian(
        self, return_order: Optional[int] = None, unit: Optional[str] = None
    ) -> Union["NDArray", "PintArray"]:
        """Get the effective Hamiltonian of the pair system.

        Args:
            return_order: The order of the perturbation to return.
                Default None, returns the sum up to the perturbation order set in the class.
            unit: The unit in which to return the effective Hamiltonian.
                If None, returns a pint array.

        Returns:
            The effective Hamiltonian of the pair system in the given unit.
            If unit is None, returns a pint array, otherwise returns a numpy array.

        """
        if self._eff_h_dict_au is None:
            self._create_effective_hamiltonian()
            assert self._eff_h_dict_au is not None
        if return_order is None:
            h_eff_au: NDArray = sum(self._eff_h_dict_au.values())  # type: ignore [assignment]
        elif return_order in self._eff_h_dict_au:
            h_eff_au = self._eff_h_dict_au[return_order]
        else:
            raise ValueError(
                f"The perturbation order {return_order} is not available in the effective Hamiltonian "
                f"with the specified perturbation_order {self.perturbation_order}."
            )
        return QuantityArray.convert_au_to_user(h_eff_au, "energy", unit)

    def get_effective_basisvectors(self) -> "csr_matrix":
        """Get the eigenvectors of the perturbative Hamiltonian."""
        if self._eff_vecs is None:
            self._create_effective_hamiltonian()
            assert self._eff_vecs is not None
        return self._eff_vecs

    def _create_effective_hamiltonian(self) -> None:
        """Calculate the perturbative Hamiltonian up to the given perturbation order."""
        hamiltonian_au = self.system_pair.get_hamiltonian(unit="hartree")
        eff_h_dict_au, eigvec_perturb = calculate_perturbative_hamiltonian(
            hamiltonian_au, self.model_inds, self.perturbation_order
        )
        self._eff_h_dict_au = eff_h_dict_au
        self._eff_vecs = eigvec_perturb

        self.check_for_resonances()

    # # # Other stuff # # #
    @cached_property
    def model_inds(self) -> list[int]:
        """The indices of the corresponding KetPairs of the given ket_tuples in the basis of the pair system."""
        model_inds = []
        for kets in self.ket_tuples:
            overlap = self.basis_pair.get_overlaps(kets)
            inds = np.argsort(overlap)[::-1]
            model_inds.append(int(inds[0]))

            if overlap[inds[0]] > 0.6:
                continue
            if overlap[inds[0]] == 0:
                raise ValueError(f"The pairstate {kets} is not part of the basis of the pair system.")
            logger.critical(
                "The ket_pair %s only has an overlap of %.3f with its corresponding pair_state."
                " The most perturbing states are:",
                kets,
                overlap[inds[0]],
            )
            for i in inds[1:5]:
                logger.error("  - %s with overlap %.3e", self.system_pair.basis.kets[i], overlap[i])

        return model_inds

    def check_for_resonances(self, required_overlap: float = 0.9) -> None:
        r"""Check if states of the model space have strong resonances with states outside the model space."""
        eigvec_perturb = self.get_effective_basisvectors()
        overlaps = (eigvec_perturb.multiply(eigvec_perturb.conj())).real  # elementwise multiplication

        model_inds = self.model_inds
        for i, m_ind in enumerate(model_inds):
            inf_data_inds = np.isinf(overlaps[i, :].data)
            if inf_data_inds.any():
                indices = overlaps[i, :].indices[np.argwhere(inf_data_inds).flatten()]
                logger.critical(
                    "Detected 'inf' entries in the effective eigenvectors.\n"
                    " This might happen, if you forgot to include a degenerate state in the model space.\n"
                    " Consider adding the following states to the model space:"
                )
                for index in indices:
                    logger.critical("  - %s has infinite admixture", self.system_pair.basis.kets[index])
                continue

            overlaps[i, :] /= sparse.linalg.norm(overlaps[i, :])
            if overlaps[i, m_ind] >= required_overlap:
                continue
            logger.error(
                "The ket %s has only %.3f overlap with its corresponding effective eigenvector.\n"
                " Thus, the calculation might lead to unexpected or wrong results.\n"
                " Consider adding the most perturbing states to the model space.\n"
                " The most perturbing states are:",
                self.system_pair.basis.kets[m_ind],
                overlaps[i, m_ind],
            )
            print_above_admixture = (1 - required_overlap) * 0.05
            indices = sparse.find(overlaps[i, :] >= print_above_admixture)[1]
            for index in indices:
                if index != m_ind:
                    admixture = overlaps[i, index]
                    logger.error("  - %s with overlap %.3e", self.system_pair.basis.kets[index], admixture)


@lru_cache(maxsize=20)
def get_basis_atom_with_cache(
    species: str, n: tuple[int, int], l: tuple[int, int], m: tuple[int, int]
) -> pi_real.BasisAtom:
    """Get a BasisAtom object potentially by using a cache to avoid recomputing it."""
    return pi_real.BasisAtom(species, n=n, l=l, m=m)
