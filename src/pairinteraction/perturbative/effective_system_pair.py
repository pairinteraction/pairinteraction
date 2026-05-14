# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import contextlib
import copy
import logging
from functools import cached_property, lru_cache
from typing import TYPE_CHECKING, Literal, overload

import numpy as np
from scipy import sparse
from typing_extensions import deprecated

from pairinteraction.basis import BasisAtom, BasisAtomReal, BasisPair, BasisPairReal
from pairinteraction.diagonalization import diagonalize
from pairinteraction.perturbative.perturbation_theory import calculate_perturbative_hamiltonian
from pairinteraction.system import SystemAtom, SystemAtomReal, SystemPair, SystemPairReal
from pairinteraction.units import QuantityArray

if TYPE_CHECKING:
    from collections.abc import Sequence

    from scipy.sparse import csr_matrix
    from typing_extensions import Self

    from pairinteraction.ket import KetAtom, KetAtomTuple  # noqa: F401
    from pairinteraction.units import ArrayLike, NDArray, PintArray, PintFloat


logger = logging.getLogger(__name__)

BasisSystemLiteral = Literal["basis_atoms", "system_atoms", "basis_pair", "system_pair"]


class EffectiveSystemPair:
    """Class for creating an effective SystemPair object and calculating the effective Hamiltonian.

    Given a subspace spanned by tuples of `KetAtom` objects (ket_tuples),
    this class automatically generates appropriate `BasisAtom`, `SystemAtom` objects as well as a `BasisPair` and
    `SystemPair` object to calculate the effective Hamiltonian in the subspace via perturbation theory.

    This class also allows to set magnetic and electric fields similar to the `SystemAtom` class,
    as well as the angle and distance between the two atoms like in the `SystemPair` class.

    Examples:
        >>> import pairinteraction as pi
        >>> ket_atoms = {
        ...     "+": pi.KetAtom("Rb", n=59, l=0, j=0.5, m=0.5),
        ...     "0": pi.KetAtom("Rb", n=58, l=1, j=1.5, m=1.5),
        ...     "-": pi.KetAtom("Rb", n=58, l=0, j=0.5, m=0.5),
        ... }
        >>> ket_tuples = [
        ...     (ket_atoms["+"], ket_atoms["-"]),
        ...     (ket_atoms["0"], ket_atoms["0"]),
        ...     (ket_atoms["-"], ket_atoms["+"]),
        ... ]
        >>> eff_system = pi.EffectiveSystemPair(ket_tuples)
        >>> eff_system = eff_system.set_distance(10, angle_degree=45, unit="micrometer")
        >>> eff_h = eff_system.get_effective_hamiltonian(unit="MHz")
        >>> eff_h -= np.eye(3) * eff_system.get_pair_energies("MHz")[1]
        >>> print(np.round(eff_h, 0), "MHz")
        [[292.   3.   0.]
         [  3.  -0.   3.]
         [  0.   3. 292.]] MHz

    """

    _basis_atom_class = BasisAtom
    _basis_pair_class = BasisPair
    _system_atom_class = SystemAtom
    _system_pair_class = SystemPair

    def __init__(self, ket_tuples: Sequence[KetAtomTuple]) -> None:
        if not all(len(ket_tuple) == 2 for ket_tuple in ket_tuples):
            raise ValueError("All ket tuples must contain exactly two kets")
        for i in range(2):
            if not all(ket_tuple[i].species == ket_tuples[0][i].species for ket_tuple in ket_tuples):
                raise ValueError(f"All kets for atom={i} must have the same species")

        # Perturbation attributes
        self._ket_tuples = [tuple(kets) for kets in ket_tuples]
        self._perturbation_order = 2

        # BasisAtom and SystemAtom attributes
        self._delta_n: int | None = None
        self._delta_l: int | None = None
        self._delta_m: int | None = None
        self._electric_field: PintArray | None = None
        self._magnetic_field: PintArray | None = None
        self._diamagnetism_enabled: bool | None = None

        # BasisPair and SystemPair attributes
        self._interaction_order: int | None = None
        self._distance_vector: PintArray | None = None

        # misc
        self._eff_h_dict_au: dict[int, NDArray] | None = None
        self._eff_vecs: csr_matrix | None = None

        # misc user set stuff
        self._user_set_parts: set[BasisSystemLiteral] = set()

    def copy(self: Self) -> Self:
        """Create a copy of the EffectiveSystemPair object (before it has been created)."""
        if self._is_created("basis_atoms"):
            raise RuntimeError(
                "Cannot copy the EffectiveSystemPair object after it has been created. "
                "Please create a new object instead."
            )
        return copy.copy(self)

    def _is_created(self: Self, what: BasisSystemLiteral = "basis_atoms") -> bool:
        """Check if some part of the effective Hamiltonian has already been created."""
        return hasattr(self, "_" + what)

    def _ensure_not_created(self: Self, what: BasisSystemLiteral = "basis_atoms") -> None:
        """Ensure that some part of the effective Hamiltonian has not been created yet."""
        if self._is_created(what):
            raise RuntimeError(
                f"Cannot change parameters for {what} after it has already been created. "
                f"Please set all parameters before {what} before accessing it (or creating the effective Hamiltonian)."
            )

    def _delete_created(self: Self, what: BasisSystemLiteral = "basis_atoms") -> None:
        """Delete the created part of the effective Hamiltonian.

        Args:
            what: The part of the effective Hamiltonian to delete.
                Default is "basis_atoms", which means delete all parts that have been created.

        """
        self._eff_h_dict_au = None
        self._eff_vecs = None
        self._eff_basis = None
        with contextlib.suppress(AttributeError):
            del self.model_inds

        parts_order: list[BasisSystemLiteral] = ["system_pair", "basis_pair", "system_atoms", "basis_atoms"]
        for part in parts_order:
            if part in self._user_set_parts:
                raise RuntimeError(
                    f"Cannot delete {part} because it has been set by the user. "
                    "Please create a new EffectiveSystemPair object instead."
                )
            with contextlib.suppress(AttributeError):
                delattr(self, "_" + part)
            if part == what:
                break

    # # # Perturbation methods and attributes # # #
    @property
    def ket_tuples(self) -> list[KetAtomTuple]:
        """The tuples of kets, which form the model space for the effective Hamiltonian."""
        return self._ket_tuples  # type: ignore [return-value]

    @property
    def perturbation_order(self) -> int:
        """The perturbation order for the effective Hamiltonian."""
        return self._perturbation_order

    def set_perturbation_order(self: Self, order: int) -> Self:
        """Set the perturbation order for the effective Hamiltonian."""
        self._delete_created()
        self._perturbation_order = order
        return self

    # # # BasisAtom methods and attributes # # #
    @property
    def basis_atoms(self) -> tuple[BasisAtom, BasisAtom]:
        """The basis objects for the single-atom systems."""
        if not self._is_created("basis_atoms"):
            self._create_basis_atoms()
        return self._basis_atoms  # type: ignore [return-value]

    @basis_atoms.setter
    def basis_atoms(self, basis_atoms: tuple[BasisAtom, BasisAtom]) -> None:
        self._ensure_not_created()
        if self._delta_n is not None or self._delta_l is not None or self._delta_m is not None:
            logger.warning("Setting basis_atoms will overwrite parameters defined for basis_atoms.")
        self._user_set_parts.add("basis_atoms")
        self._basis_atoms = tuple(basis_atoms)

    def set_delta_n(self: Self, delta_n: int) -> Self:
        """Set the delta_n value for single-atom basis."""
        self._delete_created()
        self._delta_n = delta_n
        return self

    def set_delta_l(self: Self, delta_l: int) -> Self:
        """Set the delta_l value for single-atom basis."""
        self._delete_created()
        self._delta_l = delta_l
        return self

    def set_delta_m(self: Self, delta_m: int) -> Self:
        """Set the delta_m value for single-atom basis."""
        self._delete_created()
        self._delta_m = delta_m
        return self

    def _create_basis_atoms(self) -> None:
        delta_n = self._delta_n if self._delta_n is not None else 7
        delta_l = self._delta_l
        if delta_l is None:
            delta_l = self.perturbation_order * (self.interaction_order - 2)
        delta_m = self._delta_m
        if delta_m is None and self._delta_l is None and self._are_fields_along_z:
            delta_m = self.perturbation_order * (self.interaction_order - 2)

        basis_atoms: list[BasisAtom] = []
        use_real = isinstance(self, EffectiveSystemPairReal)
        for i in range(2):
            kets = [ket_tuple[i] for ket_tuple in self.ket_tuples]
            nlfm = np.transpose([[ket.n, ket.l, ket.f, ket.m] for ket in kets])
            n_range = (int(np.min(nlfm[0])) - delta_n, int(np.max(nlfm[0])) + delta_n)
            l_range = (np.min(nlfm[1]) - delta_l, np.max(nlfm[1]) + delta_l)
            if any(ket.is_calculated_with_mqdt for ket in kets) and self._delta_l is None:
                # for mqdt we increase the default delta_l by 1 to take into account the variance ...
                l_range = (np.min(nlfm[1]) - delta_l - 1, np.max(nlfm[1]) + delta_l + 1)
            m_range = (np.min(nlfm[3]) - delta_m, np.max(nlfm[3]) + delta_m) if delta_m is not None else None
            basis = get_basis_atom_with_cache(kets[0].species, n_range, l_range, m_range, use_real=use_real)
            basis_atoms.append(basis)

        self._basis_atoms = tuple(basis_atoms)

    # # # SystemAtom methods and attributes # # #
    @property
    def system_atoms(self) -> tuple[SystemAtom, SystemAtom]:
        """The system objects for the single-atom systems."""
        if not self._is_created("system_atoms"):
            self._create_system_atoms()
        return self._system_atoms

    @system_atoms.setter
    def system_atoms(self, system_atoms: tuple[SystemAtom, SystemAtom]) -> None:
        self._ensure_not_created()
        if (
            self._electric_field is not None
            or self._magnetic_field is not None
            or self._diamagnetism_enabled is not None
        ):
            logger.warning("Setting system_atoms will overwrite parameters defined for system_atoms.")
        self._user_set_parts.add("system_atoms")
        self._system_atoms: tuple[SystemAtom, SystemAtom] = tuple(system_atoms)  # type: ignore [assignment]
        self.basis_atoms = tuple(system.basis for system in system_atoms)  # type: ignore [assignment]

    @property
    def electric_field(self) -> PintArray:
        """The electric field for the single-atom systems."""
        if self._electric_field is None:
            self.set_electric_field([0, 0, 0], "V/cm")
            assert self._electric_field is not None
        return self._electric_field

    def set_electric_field(
        self: Self,
        electric_field: PintArray | ArrayLike,
        unit: str | None = None,
    ) -> Self:
        """Set the electric field for the single-atom systems.

        Args:
            electric_field: The electric field to set for the systems.
            unit: The unit of the electric field, e.g. "V/cm".
                Default None expects a `pint.Quantity`.

        """
        self._delete_created()
        self._electric_field = QuantityArray.convert_user_to_pint(electric_field, unit, "electric_field")
        return self

    @property
    def magnetic_field(self) -> PintArray:
        """The magnetic field for the single-atom systems."""
        if self._magnetic_field is None:
            self.set_magnetic_field([0, 0, 0], "gauss")
            assert self._magnetic_field is not None
        return self._magnetic_field

    def set_magnetic_field(
        self: Self,
        magnetic_field: PintArray | ArrayLike,
        unit: str | None = None,
    ) -> Self:
        """Set the magnetic field for the single-atom systems.

        Args:
            magnetic_field: The magnetic field to set for the systems.
            unit: The unit of the magnetic field, e.g. "gauss".
                Default None expects a `pint.Quantity`.

        """
        self._delete_created()
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

    def set_diamagnetism_enabled(self: Self, enable: bool = True) -> Self:
        """Enable or disable diamagnetism for the system.

        Args:
            enable: Whether to enable or disable diamagnetism.

        """
        self._delete_created("system_atoms")
        self._diamagnetism_enabled = enable
        return self

    def _create_system_atoms(self) -> None:
        system_atoms: list[SystemAtom] = []
        for basis_atom in self.basis_atoms:
            system = self._system_atom_class(basis_atom)
            system.set_diamagnetism_enabled(self.diamagnetism_enabled)
            system.set_electric_field(self.electric_field)
            system.set_magnetic_field(self.magnetic_field)
            system_atoms.append(system)
        diagonalize(system_atoms)

        self._system_atoms = tuple(system_atoms)  # type: ignore [assignment]

    @overload
    def get_pair_energies(self, unit: None = None) -> list[PintFloat]: ...

    @overload
    def get_pair_energies(self, unit: str) -> list[float]: ...

    def get_pair_energies(self, unit: str | None = None) -> list[float] | list[PintFloat]:
        """Get the pair energies of the ket tuples for infinite distance (i.e. no interaction).

        Args:
            unit: The unit to which to convert the energies to.
                Default None will return a list of `pint.Quantity`.

        Returns:
            The energies as list of float if a unit was given, otherwise as list of `pint.Quantity`.

        """
        return [  # type: ignore [return-value]
            sum(
                system.get_corresponding_energy(ket, unit=unit)
                for system, ket in zip(self.system_atoms, ket_tuple, strict=True)
            )
            for ket_tuple in self.ket_tuples
        ]

    # # # BasisPair methods and attributes # # #
    @property
    def basis_pair(self) -> BasisPair:
        """The basis pair object for the pair system."""
        if not self._is_created("basis_pair"):
            self.create_basis_pair()
        return self._basis_pair

    @basis_pair.setter
    def basis_pair(self, basis_pair: BasisPair) -> None:
        self._ensure_not_created()
        self._user_set_parts.add("basis_pair")
        self._basis_pair = basis_pair
        self.system_atoms = basis_pair.system_atoms

    @deprecated("set_minimum_number_of_ket_pairs is deprecated, use create_basis_pair(...) instead.")
    def set_minimum_number_of_ket_pairs(self: Self, number_of_kets: int) -> Self:  # noqa: ARG002
        raise DeprecationWarning("set_minimum_number_of_ket_pairs is deprecated, use create_basis_pair(...) instead.")

    @deprecated("set_maximum_number_of_ket_pairs is deprecated, use create_basis_pair(...) instead.")
    def set_maximum_number_of_ket_pairs(self: Self, number_of_kets: int) -> Self:  # noqa: ARG002
        raise DeprecationWarning("set_maximum_number_of_ket_pairs is deprecated, use create_basis_pair(...) instead.")

    def create_basis_pair(
        self,
        delta_energy: float | PintFloat | None = None,
        delta_energy_unit: str | None = None,
        number_of_kets: int | None = None,
        *,
        allow_large_basis: bool = False,
    ) -> None:
        if self._is_created("basis_pair"):
            raise RuntimeError("The basis_pair has already been created. Cannot create it again.")

        if delta_energy is not None or number_of_kets is not None:
            self._basis_pair = self._basis_pair_class.from_kets(
                self.ket_tuples,
                system_atoms=self.system_atoms,
                delta_energy=delta_energy,
                delta_energy_unit=delta_energy_unit,
                number_of_kets=number_of_kets,
            )
            return

        min_nu = min(ket.nu for ket_tuple in self.ket_tuples for ket in ket_tuple)
        # for nu = 40 use delta_energy = 8GHz and scale with nu^3 (i.e. for nu=80 use 1GHz)
        delta_energy_ghz = 8 * (40 / min_nu) ** 3
        basis_pair = self._basis_pair_class.from_kets(
            self.ket_tuples,
            system_atoms=self.system_atoms,
            delta_energy=delta_energy_ghz,
            delta_energy_unit="GHz",
        )

        if basis_pair.number_of_kets > 25_000:
            msg = (
                f"The automatically generated basis_pair contains {basis_pair.number_of_kets} kets. "
                "This might lead to long calculation times for the effective Hamiltonian. "
            )
            if not allow_large_basis:
                raise RuntimeError(
                    msg
                    + "If this is on purpose, consider calling `create_basis_pair(allow_large_basis=True)`. "
                    + "If not, consider calling `create_basis_pair(delta_energy=..., delta_energy_unit=...)` "
                    + "or `create_basis_pair(number_of_kets=...)` "
                    + "with custom parameters to control the basis size. "
                )
            logger.warning(msg)

        self._basis_pair = basis_pair
        logger.debug("The pair basis for the perturbative calculations consists of %d kets.", basis_pair.number_of_kets)

    # # # SystemPair methods and attributes # # #
    @property
    def system_pair(self) -> SystemPair:
        """The system pair object for the pair system."""
        if not self._is_created("system_pair"):
            self._create_system_pair()
        return self._system_pair

    @system_pair.setter
    def system_pair(self, system_pair: SystemPair) -> None:
        self._ensure_not_created()
        if self._interaction_order is not None or self._distance_vector is not None:
            logger.warning("Setting system_pair will overwrite parameters defined for system_pair.")
        self._user_set_parts.add("system_pair")
        self._system_pair = system_pair
        self.basis_pair = system_pair.basis

    @property
    def interaction_order(self) -> int:
        """The interaction order for the pair system."""
        if self._interaction_order is None:
            self.set_interaction_order(3)
        return self._interaction_order  # type: ignore [return-value]

    def set_interaction_order(self: Self, order: int) -> Self:
        """Set the interaction order of the pair system.

        Args:
            order: The interaction order to set for the pair system.
                The order must be 3, 4, or 5.

        """
        self._delete_created()
        self._interaction_order = order
        return self

    @property
    def distance_vector(self) -> PintArray:
        """The distance vector between the atoms in the pair system."""
        if self._distance_vector is None:
            self.set_distance_vector([0, 0, np.inf], "micrometer")
        return self._distance_vector  # type: ignore [return-value]

    def set_distance(
        self: Self,
        distance: float | PintFloat,
        angle_degree: float = 0,
        unit: str | None = None,
    ) -> Self:
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
        self: Self,
        distance: ArrayLike | PintArray,
        unit: str | None = None,
    ) -> Self:
        """Set the distance vector between the atoms.

        Args:
            distance: The distance vector to set between the atoms in the given unit.
            unit: The unit of the distance, e.g. "micrometer".
                Default None expects a `pint.Quantity`.

        """
        self._delete_created("system_pair")
        self._distance_vector = QuantityArray.convert_user_to_pint(distance, unit, "distance")
        return self

    def set_angle(
        self: Self,
        angle: float = 0,
        unit: Literal["degree", "radian"] = "degree",
    ) -> Self:
        """Set the angle between the atoms in degrees.

        Args:
            angle: The angle between the distance vector and the z-axis (by default in degrees).
                90 degrees corresponds to the x-axis.
                Defaults to 0, which corresponds to the z-axis.
            unit: The unit of the angle, either "degree" or "radian", by default "degree".

        """
        assert unit in ("radian", "degree"), f"Unit {unit} is not supported for angle."
        if unit == "radian":
            angle = np.rad2deg(angle)
        distance_mum: float = np.linalg.norm(self.distance_vector.to("micrometer").magnitude)  # type: ignore [assignment]
        return self.set_distance(distance_mum, angle, "micrometer")

    def _create_system_pair(self) -> None:
        system_pair = self._system_pair_class(self.basis_pair)
        system_pair.set_distance_vector(self.distance_vector)
        system_pair.set_interaction_order(self.interaction_order)
        self._system_pair = system_pair

    # # # Effective Hamiltonian methods and attributes # # #
    @overload
    def get_effective_hamiltonian(self, return_order: int | None = None, unit: None = None) -> PintArray: ...

    @overload
    def get_effective_hamiltonian(self, return_order: int | None = None, *, unit: str) -> NDArray: ...

    def get_effective_hamiltonian(
        self, return_order: int | None = None, unit: str | None = None
    ) -> NDArray | PintArray:
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
        return QuantityArray.convert_au_to_user(np.real_if_close(h_eff_au), "energy", unit)

    def get_effective_basisvectors(self) -> csr_matrix:
        """Get the eigenvectors of the perturbative Hamiltonian."""
        if len(self.model_inds) > 1 and self.perturbation_order > 2:
            logger.warning("For more than one state and perturbation_order > 2 the effective basis might be wrong.")
        if self._eff_vecs is None:
            self._create_effective_hamiltonian()
            assert self._eff_vecs is not None
        return self._eff_vecs

    def get_effective_basis(self) -> BasisPair:
        """Get the effective basis of the pair system."""
        raise NotImplementedError("The get effective basis method is not implemented yet.")

    def _create_effective_hamiltonian(self) -> None:
        """Calculate the perturbative Hamiltonian up to the given perturbation order."""
        hamiltonian_au = self.system_pair.get_hamiltonian(unit="hartree")
        eff_h_dict_au, eff_vecs = calculate_perturbative_hamiltonian(
            hamiltonian_au, self.model_inds, self.perturbation_order
        )
        self._eff_h_dict_au = eff_h_dict_au
        self._eff_vecs = eff_vecs

        self.check_for_resonances()

    # # # Other stuff # # #
    @cached_property
    def model_inds(self) -> list[int]:
        """The indices of the corresponding KetPairs of the given ket_tuples in the basis_pair."""
        model_inds = []
        for kets in self.ket_tuples:
            overlap = self.basis_pair.get_overlaps(kets)
            inds = np.argsort(overlap)[::-1]
            model_inds.append(int(inds[0]))
            self._warn_model_inds_overlap(overlap, inds, kets)
        return model_inds

    def _warn_model_inds_overlap(self, overlap: NDArray, inds: NDArray, kets: KetAtomTuple) -> None:
        if overlap[inds[0]] > 0.8:
            return

        if overlap[inds[0]] == 0:
            raise ValueError(f"The pairstate {kets} is not part of the basis_pair.")

        msg = ""
        accumulated = overlap[inds[0]]
        for i in inds[1:5]:
            msg += f"\n  - {self.basis_pair.get_state(i)} with overlap {overlap[i]:.3e}"
            accumulated += overlap[i]
            if accumulated > 0.8:
                break

        logger.warning(
            "The pairstate %s has only an overlap of %.3f with its corresponding state in the basis_pair.\n"
            "Note that the effective hamiltonian is calculated with respect to the corresponding state %s.\n"
            "The most perturbing other states in the basis_pair are:\n%s",
            *(kets, overlap[inds[0]], self.basis_pair.get_state(inds[0]), msg),
        )

    def check_for_resonances(self, max_perturber_weight: float = 0.05) -> None:
        r"""Check if states of the model space have strong resonances with states outside the model space."""
        # Get the effective eigenvectors without potential warning
        if self._eff_vecs is None:
            self._create_effective_hamiltonian()
            assert self._eff_vecs is not None
        eff_vecs = self._eff_vecs

        overlaps = (eff_vecs.multiply(eff_vecs.conj())).real  # elementwise multiplication

        for i, m_ind in enumerate(self.model_inds):
            overlaps_i = overlaps[i, :]
            other_weight = np.sum(overlaps_i.data) - 1
            if other_weight < max_perturber_weight:
                continue

            msg = ""
            indices = [
                int(index) for index in sparse.find(overlaps_i >= 0.1 * max_perturber_weight)[1] if index != m_ind
            ]
            indices = sorted(indices, key=lambda index, ov=overlaps_i: ov[0, index], reverse=True)  # type: ignore [misc]
            overlap = 0
            for index in indices[:5]:
                admixture = overlaps_i[0, index]
                msg += f"\n  - {self.basis_pair.get_state(index)} has admixture {overlaps_i[0, index]:.3e}"
                overlap += admixture
                if overlap > 0.8 * other_weight:
                    break

            logger.warning(
                "The state (from the model space) %s gets a large dressing (%.3f overlap) "
                "in perturbation theory from other states from the basis_pair.\n"
                "Thus, treating these states perturbatively might not be accurate. "
                "Consider adding these states to the model space.\n"
                "The most perturbing states are:\n%s",
                *(self.basis_pair.get_state(m_ind), other_weight, msg),
            )


class EffectiveSystemPairReal(EffectiveSystemPair):
    _basis_atom_class = BasisAtomReal
    _basis_pair_class = BasisPairReal
    _system_atom_class = SystemAtomReal
    _system_pair_class = SystemPairReal


@lru_cache(maxsize=20)
def get_basis_atom_with_cache(
    species: str, n: tuple[int, int], l: tuple[int, int], m: tuple[int, int], *, use_real: bool
) -> BasisAtom:
    """Get a BasisAtom object potentially by using a cache to avoid recomputing it."""
    if use_real:
        return BasisAtomReal(species, n=n, l=l, m=m)
    return BasisAtom(species, n=n, l=l, m=m)
