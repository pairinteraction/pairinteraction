# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from contextlib import contextmanager
from functools import cached_property
from typing import TYPE_CHECKING, Literal, overload

import numpy as np
from scipy import sparse

from pairinteraction.basis import BasisAtom, BasisAtomReal, BasisPair, BasisPairReal
from pairinteraction.diagonalization import diagonalize
from pairinteraction.perturbative.perturbation_theory import calculate_perturbative_hamiltonian
from pairinteraction.system import SystemAtom, SystemAtomReal, SystemPair, SystemPairReal
from pairinteraction.units import QuantityArray

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    from scipy.sparse import csr_matrix
    from typing_extensions import Self

    from pairinteraction.ket import KetAtomTuple
    from pairinteraction.units import ArrayLike, NDArray, PintArray, PintFloat


logger = logging.getLogger(__name__)


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
         [  3.   0.   3.]
         [  0.   3. 292.]] MHz

    """

    _basis_atom_class = BasisAtom
    _basis_pair_class = BasisPair
    _system_atom_class = SystemAtom
    _system_pair_class = SystemPair

    def __init__(self, ket_tuples: Sequence[KetAtomTuple], perturbation_order: int = 2) -> None:
        if not all(len(ket_tuple) == 2 for ket_tuple in ket_tuples):
            raise ValueError("All ket tuples must contain exactly two kets")
        for i in range(2):
            if not all(ket_tuple[i].species == ket_tuples[0][i].species for ket_tuple in ket_tuples):
                raise ValueError(f"All kets for atom={i} must have the same species")

        # Perturbation attributes
        self._ket_tuples: list[KetAtomTuple] = [tuple(kets) for kets in ket_tuples]
        self._perturbation_order = perturbation_order

        # Basis and System attributes
        self._basis_atoms: tuple[BasisAtom, BasisAtom] | None = None
        self._system_atoms: tuple[SystemAtom, SystemAtom] | None = None
        self._basis_pair: BasisPair | None = None
        self._system_pair: SystemPair | None = None

        # System parameters
        self._electric_field: PintArray | None = None
        self._magnetic_field: PintArray | None = None
        self._diamagnetism_enabled: bool | None = None
        self._distance_vector: PintArray | None = None
        self._interaction_order: int | None = None

        # effective Hamiltonian attributes
        self._eff_h_dict_au: dict[int, NDArray] | None = None
        self._eff_vecs: csr_matrix | None = None

    @property
    def ket_tuples(self) -> list[KetAtomTuple]:
        """The tuples of kets, which form the model space for the effective Hamiltonian."""
        return self._ket_tuples

    @property
    def perturbation_order(self) -> int:
        """The perturbation order for the effective Hamiltonian."""
        return self._perturbation_order

    @property
    def basis_atoms(self) -> tuple[BasisAtom, BasisAtom]:
        """The basis objects for the single-atom systems."""
        if self._basis_atoms is None:
            self._auto_create_basis_atoms()
            assert self._basis_atoms is not None
        return self._basis_atoms

    @property
    def system_atoms(self) -> tuple[SystemAtom, SystemAtom]:
        """The system objects for the single-atom systems."""
        if self._system_atoms is None:
            self._auto_create_system_atoms()
            assert self._system_atoms is not None
        return self._system_atoms

    @property
    def basis_pair(self) -> BasisPair:
        """The basis pair object for the pair system."""
        if self._basis_pair is None:
            self.create_basis_pair()
            assert self._basis_pair is not None
        return self._basis_pair

    @property
    def system_pair(self) -> SystemPair:
        """The system pair object for the pair system."""
        if self._system_pair is None:
            self._auto_create_system_pair()
            assert self._system_pair is not None
        return self._system_pair

    def _ensure_basis_atoms_not_created(self, name: str) -> None:
        if self._basis_atoms is not None:
            raise RuntimeError(
                f"Cannot set {name} because some basis or system has already been set or created. "
                f"If you want to set {name} manually, do this directly after creating the EffectiveSystemPair object."
            )

    def _ensure_no_system_atom_parameters_are_set(self) -> None:
        if not (self._electric_field is None and self._magnetic_field is None and self._diamagnetism_enabled is None):
            raise RuntimeError("Cannot set system_atoms because some system atom parameters have already been set.")

    def _ensure_no_system_pair_parameters_are_set(self) -> None:
        if not (self._distance_vector is None and self._interaction_order is None):
            raise RuntimeError("Cannot set system_pair because some system pair parameters have already been set.")

    def set_basis_atoms(self, basis_atoms: Sequence[BasisAtom] | BasisAtom) -> None:
        self._ensure_basis_atoms_not_created("basis_atoms")
        if isinstance(basis_atoms, BasisAtom):
            basis_atoms = (basis_atoms, basis_atoms)
        if not len(basis_atoms) == 2 or not all(isinstance(basis_atom, BasisAtom) for basis_atom in basis_atoms):
            raise ValueError("basis_atoms must be a tuple of two BasisAtom objects or a single BasisAtom object.")
        self._basis_atoms = tuple(basis_atoms)  # type: ignore [assignment]

    def set_system_atoms(self, system_atoms: Sequence[SystemAtom] | SystemAtom) -> None:
        self._ensure_basis_atoms_not_created("system_atoms")
        self._ensure_no_system_atom_parameters_are_set()
        if isinstance(system_atoms, SystemAtom):
            system_atoms = (system_atoms, system_atoms)
        if not len(system_atoms) == 2 or not all(isinstance(system_atom, SystemAtom) for system_atom in system_atoms):
            raise ValueError("system_atoms must be a tuple of two SystemAtom objects or a single SystemAtom object.")
        self._system_atoms = tuple(system_atoms)  # type: ignore [assignment]
        self.set_basis_atoms([system.basis for system in system_atoms])

    def set_basis_pair(self, basis_pair: BasisPair) -> None:
        self._ensure_basis_atoms_not_created("basis_pair")
        self._basis_pair = basis_pair
        self.set_system_atoms(basis_pair.system_atoms)

    def set_system_pair(self, system_pair: SystemPair) -> None:
        self._ensure_basis_atoms_not_created("system_pair")
        self._ensure_no_system_pair_parameters_are_set()
        self._system_pair = system_pair
        self.set_basis_pair(system_pair.basis)

    def _auto_create_basis_atoms(self) -> None:
        delta_n = 7
        interaction_order = 3 if self._interaction_order is None else self._interaction_order
        delta_l = self.perturbation_order * (interaction_order - 2)

        fields_along_z = True
        if self._magnetic_field is not None:
            fields_along_z = fields_along_z and all(x == 0 for x in self._magnetic_field.m[:2])
        if self._electric_field is not None:
            fields_along_z = fields_along_z and all(x == 0 for x in self._electric_field.m[:2])

        delta_m = None
        if fields_along_z:
            delta_m = self.perturbation_order * (interaction_order - 2)

        basis_atoms: list[BasisAtom] = []
        for i in range(2):
            kets = [ket_tuple[i] for ket_tuple in self.ket_tuples]
            # for mqdt we increase the default delta_l by 1 to take into account the variance ...
            _delta_l = delta_l + 1 if any(ket.is_calculated_with_mqdt for ket in kets) else delta_l
            basis_atoms.append(
                self._basis_atom_class.from_kets(kets, delta_n=delta_n, delta_l=_delta_l, delta_m=delta_m)
            )
        self._basis_atoms = tuple(basis_atoms)  # type: ignore [assignment]

    def _auto_create_system_atoms(self) -> None:
        system_atoms: list[SystemAtom] = []
        for basis_atom in self.basis_atoms:
            system = self._system_atom_class(basis_atom)
            if self._diamagnetism_enabled is not None:
                system.set_diamagnetism_enabled(self._diamagnetism_enabled)
            if self._electric_field is not None:
                system.set_electric_field(self._electric_field)
            if self._magnetic_field is not None:
                system.set_magnetic_field(self._magnetic_field)
            system_atoms.append(system)
        diagonalize(system_atoms)

        self._system_atoms = tuple(system_atoms)  # type: ignore [assignment]

    def create_basis_pair(
        self,
        delta_energy: float | PintFloat | None = None,
        delta_energy_unit: str | None = None,
        number_of_kets: int | None = None,
        *,
        allow_large_basis: bool = False,
    ) -> None:
        if any(param is not None for param in (delta_energy, delta_energy_unit, number_of_kets)):
            basis_pair = self._basis_pair_class.from_ket_atoms(
                self.ket_tuples,
                system_atoms=self.system_atoms,
                delta_energy=delta_energy,
                delta_energy_unit=delta_energy_unit,
                number_of_kets=number_of_kets,
            )
        else:
            min_n = min(ket.n for ket_tuple in self.ket_tuples for ket in ket_tuple)
            # for n = 40 use delta_energy = 1GHz and scale with n^3
            delta_energy_ghz = 1.0 * (40) ** 3 / min_n**3
            basis_pair = self._basis_pair_class.from_ket_atoms(
                self.ket_tuples,
                system_atoms=self.system_atoms,
                delta_energy=delta_energy_ghz,
                delta_energy_unit="GHz",
            )

            min_number_of_kets = 5_000
            if basis_pair.number_of_kets < min_number_of_kets:
                with suppress_logger_warnings("pairinteraction.basis.basis_pair"):
                    basis_pair = self._basis_pair_class.from_ket_atoms(
                        self.ket_tuples,
                        system_atoms=self.system_atoms,
                        number_of_kets=min_number_of_kets,
                    )

        if basis_pair.number_of_kets > 20_000:
            msg = (
                f"The (automatically) generated basis_pair contains {basis_pair.number_of_kets} kets, "
                "which might lead to long calculation times for the effective Hamiltonian."
            )
            if not allow_large_basis:
                raise RuntimeError(
                    msg
                    + " If this is on purpose, consider calling `create_basis_pair(allow_large_basis=True)`"
                    + " If not, consider calling create_basis_pair with custom parameters to control the basis size."
                )
            logger.warning(msg)

        self._basis_pair = basis_pair
        logger.debug("The pair basis for the perturbative calculations consists of %d kets.", basis_pair.number_of_kets)

    def _auto_create_system_pair(self) -> None:
        system_pair = self._system_pair_class(self.basis_pair)
        if self._distance_vector is not None:
            system_pair.set_distance_vector(self._distance_vector)
        if self._interaction_order is not None:
            system_pair.set_interaction_order(self._interaction_order)
        self._system_pair = system_pair

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
        self._ensure_basis_atoms_not_created("electric_field")
        self._electric_field = QuantityArray.convert_user_to_pint(electric_field, unit, "electric_field")
        return self

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
        self._ensure_basis_atoms_not_created("magnetic_field")
        self._magnetic_field = QuantityArray.convert_user_to_pint(magnetic_field, unit, "magnetic_field")
        return self

    def set_diamagnetism_enabled(self: Self, enable: bool = True) -> Self:
        """Enable or disable diamagnetism for the system.

        Args:
            enable: Whether to enable or disable diamagnetism.

        """
        self._ensure_basis_atoms_not_created("diamagnetism_enabled")
        self._diamagnetism_enabled = enable
        return self

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

    def set_angle(
        self: Self,
        angle: float,
        unit: Literal["degree", "radian"] = "degree",
    ) -> Self:
        """Set the angle between the atoms in degrees.

        Args:
            angle: The angle between the distance vector and the z-axis (by default in degrees).
                90 degrees corresponds to the x-axis.
            unit: The unit of the angle, either "degree" or "radian", by default "degree".

        """
        if unit not in ("radian", "degree"):
            raise ValueError(f"Unit {unit} is not supported for angle.")
        if unit == "radian":
            angle = np.rad2deg(angle)
        if self._distance_vector is None:
            raise ValueError("You first need to set a distance before setting the angle.")
        distance_mum = np.linalg.norm(self._distance_vector.to("micrometer").magnitude)
        return self.set_distance(float(distance_mum), angle, "micrometer")

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
        self._ensure_basis_atoms_not_created("distance_vector")
        self._distance_vector = QuantityArray.convert_user_to_pint(distance, unit, "distance")
        return self

    def set_interaction_order(self: Self, order: int) -> Self:
        """Set the interaction order of the pair system.

        Args:
            order: The interaction order to set for the pair system.
                The order must be 3, 4, or 5.

        """
        self._ensure_basis_atoms_not_created("interaction_order")
        self._interaction_order = order
        return self

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
            self._calculate_effective_hamiltonian()
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
            self._calculate_effective_hamiltonian()
            assert self._eff_vecs is not None
        return self._eff_vecs

    def get_effective_basis(self) -> BasisPair:
        """Get the effective basis of the pair system."""
        raise NotImplementedError("The get effective basis method is not implemented yet.")

    def _calculate_effective_hamiltonian(self) -> None:
        """Calculate the perturbative Hamiltonian up to the given perturbation order."""
        hamiltonian_au = self.system_pair.get_hamiltonian(unit="hartree")
        eff_h_dict_au, eff_vecs = calculate_perturbative_hamiltonian(
            hamiltonian_au, self.model_inds, self.perturbation_order
        )
        self._eff_h_dict_au = eff_h_dict_au
        self._eff_vecs = eff_vecs

        self.check_for_resonances()

    @cached_property
    def model_inds(self) -> list[int]:
        """The indices of the corresponding KetPairs of the given ket_tuples in the basis_pair."""
        model_inds = []
        for kets in self.ket_tuples:
            overlap = self.basis_pair.get_overlaps(kets)
            inds = np.argsort(overlap)[::-1]
            model_inds.append(int(inds[0]))
            self._warn_model_inds(overlap, inds, kets)
        return model_inds

    def _warn_model_inds(self, overlap: NDArray, inds: NDArray, kets: KetAtomTuple) -> None:
        if overlap[inds[0]] > 0.8:
            return

        if overlap[inds[0]] == 0:
            raise ValueError(f"The pairstate {kets} is not part of the basis_pair.")

        msg = ""
        accumulated = overlap[inds[0]]
        for i in inds[1:5]:
            msg += f"\n  - {self.basis_pair.kets[i]} with overlap {overlap[i]:.3e}"
            accumulated += overlap[i]
            if accumulated > 0.8:
                break

        logger.warning(
            "The pairstate %s has only an overlap of %.3f with its corresponding state in the basis_pair.\n"
            "Note that the effective hamiltonian is calculated with respect to the corresponding state %s.\n"
            "The most perturbing states are:\n%s",
            *(kets, overlap[inds[0]], self.basis_pair.kets[inds[0]], msg),
        )

    def check_for_resonances(self, max_perturber_weight: float = 0.05) -> None:
        r"""Check if states of the model space have strong resonances with states outside the model space."""
        # Get the effective eigenvectors without potential warning
        if self._eff_vecs is None:
            self._calculate_effective_hamiltonian()
            assert self._eff_vecs is not None
        eff_vecs = self._eff_vecs

        overlaps = (eff_vecs.multiply(eff_vecs.conj())).real  # elementwise multiplication

        model_inds = self.model_inds
        for i, m_ind in enumerate(model_inds):
            overlaps_i = overlaps[i, :]
            other_weight = np.sum(overlaps_i.data) - 1
            if other_weight < max_perturber_weight:
                continue

            msg = ""
            indices = [index for index in sparse.find(overlaps_i >= 0.1 * max_perturber_weight)[1] if index != m_ind]
            indices = sorted(indices, key=lambda index, ov=overlaps_i: ov[0, index], reverse=True)  # type: ignore [misc]
            overlap = 0
            for index in indices[:5]:
                admixture = overlaps_i[0, index]
                msg += f"\n  - {self.basis_pair.kets[index]} has admixture {overlaps_i[0, index]:.3e}"
                overlap += admixture
                if overlap > 0.6 * other_weight:
                    break

            logger.warning(
                "The ket %s gets a large dressing (%.3f overlap) from other states from the basis_pair.\n"
                "Thus, treating these states perturbatively might not be accurate. "
                "Consider adding these states to the model space.\n"
                "The most perturbing states are:\n%s",
                *(self.basis_pair.kets[m_ind], other_weight, msg),
            )


class EffectiveSystemPairReal(EffectiveSystemPair):
    _basis_atom_class = BasisAtomReal
    _basis_pair_class = BasisPairReal
    _system_atom_class = SystemAtomReal
    _system_pair_class = SystemPairReal


@contextmanager
def suppress_logger_warnings(log: logging.Logger | str) -> Iterator[None]:
    if isinstance(log, str):
        log = logging.getLogger(log)
    old_level = log.level
    log.setLevel(logging.ERROR)
    try:
        yield
    finally:
        log.setLevel(old_level)
