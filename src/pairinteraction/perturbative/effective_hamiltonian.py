# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from collections.abc import Sequence
from functools import cached_property
from typing import TYPE_CHECKING, Any, Optional, Union, overload

import numpy as np

import pairinteraction.complex as pi_complex
import pairinteraction.real as pi_real
from pairinteraction.perturbative.perturbation_theory import calculate_perturbative_hamiltonian
from pairinteraction.units import AtomicUnits, QuantityArray, QuantityScalar, ureg

if TYPE_CHECKING:
    import numpy.typing as npt
    from typing_extensions import Self

    from pairinteraction._wrapped.basis.basis_atom import BasisAtom
    from pairinteraction._wrapped.basis.basis_pair import BasisPair
    from pairinteraction._wrapped.ket.ket_pair import KetAtomTuple
    from pairinteraction._wrapped.system.system_atom import SystemAtom
    from pairinteraction._wrapped.system.system_pair import SystemPair
    from pairinteraction.units import ArrayLike, NDArray, PintArray, PintFloat


logger = logging.getLogger(__name__)


class EffectiveHamiltonian:
    def __init__(self, ket_tuples: Sequence["KetAtomTuple"]) -> None:
        if not all(len(ket_tuple) == 2 for ket_tuple in ket_tuples):
            raise ValueError("All ket tuples must contain exactly two kets")
        for i in range(2):
            if not all(ket_tuple[i].species == ket_tuples[0][i].species for ket_tuple in ket_tuples):
                raise ValueError(f"All kets for atom={i} must have the same species")

        self.ket_tuples = ket_tuples
        self.effective_dim = len(ket_tuples)

        self.perturbation_order = 2
        self._started_creating = False

        self.delta_n = 7
        self._delta_l: Optional[int] = None

        self.electric_field: PintArray = ureg.Quantity([0, 0, 0], "V/cm")
        self.magnetic_field: PintArray = ureg.Quantity([0, 0, 0], "G")
        self.diamagnetism_enabled: bool = False

        self.interaction_order = 3
        self.distance_vector: PintArray = ureg.Quantity([0, 0, np.inf], "micrometer")

        self._basis_atoms: Optional[tuple[BasisAtom[Any], BasisAtom[Any]]] = None
        self._system_atoms: Optional[tuple[SystemAtom[Any], SystemAtom[Any]]] = None
        self._basis_pair: Optional[BasisPair[Any, Any]] = None
        self._system_pair: Optional[SystemPair[Any]] = None
        self._h_eff_dict: Optional[dict[int, NDArray]] = None

    @property
    def pi(self) -> Any:
        return pi_real if self.is_real else pi_complex

    @property
    def is_real(self) -> bool:
        return self.magnetic_field[1] == 0 and self.electric_field[1] == 0  # type: ignore [index]

    @property
    def are_fields_along_z(self) -> bool:
        return all(x == 0 for x in [*self.magnetic_field[:2], *self.electric_field[:2]])  # type: ignore [index]

    @property
    def delta_l(self) -> int:
        if self._delta_l is None:
            self._delta_l = self.perturbation_order * (self.interaction_order - 2)
        return self._delta_l

    @property
    def pair_energies(self) -> list["PintFloat"]:
        return [  # type: ignore [return-value]
            sum(system.get_corresponding_energy(ket) for system, ket in zip(self.system_atoms, ket_tuple))  # type: ignore [misc]
            for ket_tuple in self.ket_tuples
        ]

    @property
    def basis_atoms(self) -> tuple["BasisAtom[Any]", "BasisAtom[Any]"]:
        if self._basis_atoms is None:
            self.create_basis_atoms()
        return self._basis_atoms  # type: ignore [return-value]

    @property
    def system_atoms(self) -> tuple["SystemAtom[Any]", "SystemAtom[Any]"]:
        if self._system_atoms is None:
            self.create_system_atoms()
        return self._system_atoms  # type: ignore [return-value]

    @property
    def basis_pair(self) -> "BasisPair[Any, Any]":
        if self._basis_pair is None:
            self.create_basis_pair_from_number_of_kets()
        return self._basis_pair  # type: ignore [return-value]

    @property
    def system_pair(self) -> "SystemPair[Any]":
        if self._system_pair is None:
            self.create_system_pair()
        return self._system_pair  # type: ignore [return-value]

    @cached_property
    def model_inds(self) -> "npt.NDArray[np.int64]":
        """The indices of the corresponding KetPairs of the given ket_tuples in the basis of the pair system."""
        model_inds = []
        for kets in self.ket_tuples:
            overlap = self.basis_pair.get_overlaps(kets)
            index = np.argmax(overlap)
            if overlap[index] == 0:
                raise ValueError(f"The pairstate {kets} is not part of the basis of the pair system.")
            if overlap[index] < 0.5:
                raise ValueError(f"The pairstate {kets} cannot be identified uniquely (max overlap: {overlap[index]}).")
            model_inds.append(int(index))
        return np.array(model_inds)

    @cached_property
    def other_inds(self) -> "npt.NDArray[np.int64]":
        return np.setdiff1d(np.arange(self.basis_pair.number_of_states), self.model_inds)

    def set_basis_atoms(self, basis_atoms: tuple["BasisAtom[Any]", "BasisAtom[Any]"]) -> None:
        self._started_creating = True
        self._basis_atoms = basis_atoms

    def create_basis_atoms(self) -> None:
        self._started_creating = True

        basis_atoms: list[BasisAtom[Any]] = []
        for i in range(2):
            kets = [ket_tuple[i] for ket_tuple in self.ket_tuples]
            nlfm = np.transpose([[ket.n, ket.l, ket.f, ket.m] for ket in kets])
            n_range = (int(np.min(nlfm[0])) - self.delta_n, int(np.max(nlfm[0])) + self.delta_n)
            l_range = (np.min(nlfm[1]) - self.delta_l, np.max(nlfm[1]) + self.delta_l)
            if any(ket.is_calculated_with_mqdt for ket in kets):
                # for mqdt we increase self.delta_l by 1 to take into account the variance ...
                l_range = (np.min(nlfm[1]) - self.delta_l - 1, np.max(nlfm[1]) + self.delta_l + 1)
            m_range = None
            if self.are_fields_along_z:
                m_range = (np.min(nlfm[3]) - self.delta_l, np.max(nlfm[3]) + self.delta_l)
            basis = pi_real.BasisAtom(kets[0].species, n=n_range, l=l_range, m=m_range)
            basis_atoms.append(basis)

        self._basis_atoms = (basis_atoms[0], basis_atoms[1])

    def set_system_atoms(self, system_atoms: tuple["SystemAtom[Any]", "SystemAtom[Any]"]) -> None:
        self._started_creating = True
        self._system_atoms = system_atoms
        self._basis_atoms = tuple(system.basis for system in system_atoms)

    def create_system_atoms(self) -> None:
        self._started_creating = True

        system_atoms: list[SystemAtom[Any]] = []
        for basis_atom in self.basis_atoms:
            system = pi_real.SystemAtom(basis_atom)  # type: ignore [arg-type]
            system.set_diamagnetism_enabled(self.diamagnetism_enabled)
            system.set_electric_field(self.electric_field)
            system.set_magnetic_field(self.magnetic_field)
            system_atoms.append(system)
        pi_real.diagonalize(system_atoms)

        self._system_atoms = (system_atoms[0], system_atoms[1])

    def create_basis_pair_from_number_of_kets(self, number_of_kets: int = 2_000, *, exact: bool = False) -> None:
        pair_energies_au = [energy.to_base_units().magnitude for energy in self.pair_energies]

        def get_basis_pair(delta_energy_au: float) -> "BasisPair[Any, Any]":
            return pi_real.BasisPair(
                self.system_atoms,
                energy=(min(pair_energies_au) - delta_energy_au, max(pair_energies_au) + delta_energy_au),
                energy_unit=str(AtomicUnits["energy"]),
            )

        mhz_au = QuantityScalar.convert_user_to_au(1, "MHz", "energy")
        delta_energy_au = 100 * mhz_au
        min_delta, max_delta = None, None

        # make a bisect search to get a sensible basis size between:
        # number_of_kets and 1.2 * number_of_kets
        basis_pair = get_basis_pair(delta_energy_au)
        while 0.1 * mhz_au < delta_energy_au < 1:
            if max_delta is not None and min_delta is not None and max_delta - min_delta < mhz_au:
                break  # stop condition if delta_energy_au does not change anymore

            if basis_pair.number_of_kets < number_of_kets:
                min_delta = delta_energy_au
                if max_delta is None:
                    delta_energy_au *= 2
                else:
                    delta_energy_au = (delta_energy_au + max_delta) / 2
            elif basis_pair.number_of_kets > number_of_kets * (1 if exact else 1.2):
                max_delta = delta_energy_au
                if min_delta is None:
                    delta_energy_au /= 2
                else:
                    delta_energy_au = (delta_energy_au + min_delta) / 2
            else:
                break
            basis_pair = get_basis_pair(delta_energy_au)

        self._basis_pair = basis_pair
        logger.debug("The pair basis for the perturbative calculations consists of %d kets.", basis_pair.number_of_kets)

    def create_system_pair(self) -> None:
        self._system_pair = pi_real.SystemPair(self.basis_pair)  # type: ignore [arg-type]
        self._system_pair.set_distance_vector(self.distance_vector)
        self._system_pair.set_interaction_order(self.interaction_order)

    def create_effective_hamiltonian(self) -> None:
        """Calculate the perturbative Hamiltonian up to the given perturbation order."""
        hamiltonian_au = self.system_pair.get_hamiltonian(unit=str(AtomicUnits["energy"]))
        h_eff_dict, eigvec_perturb = calculate_perturbative_hamiltonian(
            hamiltonian_au, self.model_inds, self.other_inds, self.perturbation_order
        )
        self._h_eff_dict = h_eff_dict
        self._eigvec_perturb = eigvec_perturb

    def set_perturbation_order(self: "Self", order: int) -> "Self":
        self.perturbation_order = order
        return self

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
        self.electric_field = QuantityArray.from_pint_or_unit(electric_field, unit, "electric_field").to_pint()
        return self

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
        self.magnetic_field = QuantityArray.from_pint_or_unit(magnetic_field, unit, "magnetic_field").to_pint()
        return self

    def set_diamagnetism_enabled(self: "Self", enable: bool = True) -> "Self":
        self.diamagnetism_enabled = enable
        return self

    def set_interaction_order(self: "Self", order: int) -> "Self":
        self.order = order
        return self

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
        self.distance_vector = QuantityArray.from_pint_or_unit(distance, unit, "distance").to_pint()
        return self

    @overload
    def get_effective_hamiltonian(self, unit: None = None) -> "PintArray": ...

    @overload
    def get_effective_hamiltonian(self, unit: str) -> "NDArray": ...

    def get_effective_hamiltonian(self, unit: Optional[str] = None) -> Union["NDArray", "PintArray"]:
        if self._h_eff_dict is None:
            self.create_effective_hamiltonian()
            assert self._h_eff_dict is not None
        h_eff_au: NDArray = sum(h for h in self._h_eff_dict.values())  # type: ignore [assignment]
        return QuantityArray.convert_au_to_user(h_eff_au, "energy", unit)

    def create_resonances_lists(self) -> None:
        # FIXME: it is rather slow to build the full Hamiltonian
        # optimize this at some point to only calculate the necessary rows
        self.system_pair.get_hamiltonian()
        hamiltonian_au = self.system_pair.get_hamiltonian("hartree")
        energies_au = hamiltonian_au.diagonal()

        ket_pairs = self.system_pair.basis.kets
        gaps_au = np.array([energies_au - pair_energy.to_base_units().magnitude for pair_energy in self.pair_energies])

        # set all elements in the model space to zero, to only get the couplings to the other space
        hamiltonian_au[np.ix_(self.model_inds, self.model_inds)] = 0
        couplings_au = np.array([hamiltonian_au.getrow(m_ind).toarray()[0] for m_ind in self.model_inds])

        ev_contributions = np.abs(couplings_au / gaps_au)
        ev_contributions_max = np.max(ev_contributions, axis=0)

        inds = np.argsort(ev_contributions_max)[::-1]
        inds = np.array([i for i in inds if ev_contributions_max[i] > 0])

        self._all_ket_pairs = np.array(ket_pairs)[inds].tolist()
        self._all_gaps_au = gaps_au[:, inds]
        self._all_couplings_au = couplings_au[:, inds]
        self._all_ev_contributions_max = ev_contributions_max[inds]
