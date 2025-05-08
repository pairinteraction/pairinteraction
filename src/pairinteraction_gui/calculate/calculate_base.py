# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from abc import ABC
from typing import TYPE_CHECKING, Any, Generic, Optional, TypeVar, Union

import numpy as np
from attr import dataclass

from pairinteraction import (
    _wrapped,
    complex as pi_complex,
    real as pi_real,
)
from pairinteraction_gui.config.system_config import RangesKeys

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.units import NDArray
    from pairinteraction_gui.page import OneAtomPage, TwoAtomsPage

logger = logging.getLogger(__name__)

# FIXME: having all kwargs dictionaries being Any is a hacky solution, it would be nice to use TypedDict in the future

UnitFromRangeKey: dict[RangesKeys, str] = {
    "Ex": "V/cm",
    "Ey": "V/cm",
    "Ez": "V/cm",
    "Bx": "Gauss",
    "By": "Gauss",
    "Bz": "Gauss",
    "Distance": r"$\mu$m",
    "Angle": r"$^\circ$",
}

VariableNameFromRangeKey: dict[RangesKeys, str] = {
    "Ex": "efield_x",
    "Ey": "efield_y",
    "Ez": "efield_z",
    "Bx": "bfield_x",
    "By": "bfield_y",
    "Bz": "bfield_z",
    "Distance": "distance",
    "Angle": "angle",
}

PageType = TypeVar("PageType", "OneAtomPage", "TwoAtomsPage")


@dataclass
class Parameters(ABC, Generic[PageType]):
    species: tuple[str, ...]
    quantum_numbers: tuple[dict[str, float], ...]
    quantum_number_deltas: tuple[dict[str, float], ...]
    ranges: dict[RangesKeys, list[float]]
    diagonalize_kwargs: dict[str, str]
    diagonalize_relative_energy_range: Union[tuple[float, float], None]
    number_state_labels: int

    def __post_init__(self) -> None:
        """Post-initialization processing."""
        # Check if all ranges have the same number of steps
        if not all(len(v) == self.steps for v in self.ranges.values()):
            raise ValueError("All ranges must have the same number of steps")

        # Check if all tuples have the same length
        if not all(
            len(tup) == self.n_atoms for tup in [self.species, self.quantum_numbers, self.quantum_number_deltas]
        ):
            raise ValueError("All tuples must have the same length as the number of atoms")

    @classmethod
    def from_page(cls, page: PageType) -> "Self":
        """Create Parameters object from page."""
        n_atoms = page.ket_config.n_atoms

        species = tuple(page.ket_config.get_species(atom) for atom in range(n_atoms))
        quantum_numbers = tuple(page.ket_config.get_quantum_numbers(atom) for atom in range(n_atoms))

        quantum_number_deltas = tuple(page.basis_config.get_quantum_number_deltas(atom) for atom in range(n_atoms))

        ranges = page.system_config.get_ranges_dict()

        diagonalize_kwargs = {}
        if page.calculation_config.fast_mode.isChecked():
            diagonalize_kwargs["diagonalizer"] = "lapacke_evr"
            diagonalize_kwargs["float_type"] = "float32"

        diagonalize_relative_energy_range = None
        if page.calculation_config.energy_range.isChecked():
            diagonalize_relative_energy_range = page.calculation_config.energy_range.values()

        return cls(
            species,
            quantum_numbers,
            quantum_number_deltas,
            ranges,
            diagonalize_kwargs,
            diagonalize_relative_energy_range,
            page.calculation_config.number_state_labels.value(default=0),
        )

    @property
    def is_real(self) -> bool:
        """Check if the parameters are real."""
        return all(e == 0 for e in self.ranges.get("Ey", [0])) and all(b == 0 for b in self.ranges.get("By", [0]))

    @property
    def steps(self) -> int:
        """Return the number of steps."""
        return len(next(iter(self.ranges.values())))

    @property
    def n_atoms(self) -> int:
        """Return the number of atoms."""
        return len(self.species)

    def get_efield(self, step: int) -> list[float]:
        """Return the electric field for the given step."""
        efield_keys: list[RangesKeys] = ["Ex", "Ey", "Ez"]
        return [self.ranges[key][step] if key in self.ranges else 0 for key in efield_keys]

    def get_bfield(self, step: int) -> list[float]:
        """Return the magnetic field for the given step."""
        bfield_keys: list[RangesKeys] = ["Bx", "By", "Bz"]
        return [self.ranges[key][step] if key in self.ranges else 0 for key in bfield_keys]

    def get_species(self, atom: Optional[int] = None) -> str:
        """Return the species for the given ket."""
        return self.species[self._check_atom(atom)]

    def get_quantum_numbers(self, atom: Optional[int] = None) -> dict[str, Any]:
        """Return the quantum numbers for the given ket."""
        return self.quantum_numbers[self._check_atom(atom)]

    def get_quantum_number_restrictions(self, atom: Optional[int] = None) -> dict[str, Any]:
        """Return the quantum number restrictions for the given ket."""
        atom = self._check_atom(atom)
        qn_restrictions: dict[str, tuple[float, float]] = {}
        for key, delta in self.quantum_number_deltas[atom].items():
            if key in self.quantum_numbers[atom]:
                qn_restrictions[key] = (
                    self.quantum_numbers[atom][key] - delta,
                    self.quantum_numbers[atom][key] + delta,
                )
            else:
                raise ValueError(f"Quantum number delta {key} not found in quantum numbers.")
        return qn_restrictions

    def _check_atom(self, atom: Optional[int] = None) -> int:
        """Check if the atom is valid."""
        if atom is not None:
            return atom
        if self.n_atoms == 1:
            return 0
        raise ValueError("Atom index is required for multiple atoms")

    def get_diagonalize_energy_range(self, energy_of_interest: float) -> dict[str, Any]:
        """Return the kwargs for the diagonalization energy range."""
        if self.diagonalize_relative_energy_range is None:
            return {}
        kwargs: dict[str, Any] = {"energy_unit": "GHz"}
        kwargs["energy_range"] = (
            energy_of_interest + self.diagonalize_relative_energy_range[0],
            energy_of_interest + self.diagonalize_relative_energy_range[1],
        )
        return kwargs

    def get_x_values(self) -> list[float]:
        """Return the x values for the plot."""
        max_key = self._get_ranges_max_diff_key()
        return self.ranges[max_key]

    def get_x_label(self) -> str:
        """Return the x values for the plot."""
        max_key = self._get_ranges_max_diff_key()
        x_label = f"{max_key} [{UnitFromRangeKey[max_key]}]"

        non_constant_keys = [key for key, values in self.ranges.items() if key != max_key and values[0] != values[-1]]
        if non_constant_keys:
            x_label += f"  ({', '.join(non_constant_keys)} did also change)"

        return x_label

    def _get_ranges_max_diff_key(self) -> RangesKeys:
        """Return the key with the maximum difference in the ranges."""
        range_diffs: dict[RangesKeys, float] = {key: abs(r[-1] - r[0]) for key, r in self.ranges.items()}
        return max(range_diffs, key=lambda x: range_diffs.get(x, -1))

    def to_replacement_dict(self) -> dict[str, str]:
        """Return a dictionary with the parameters for replacement."""
        max_key = self._get_ranges_max_diff_key()
        replacements: dict[str, str] = {
            "$PI_DTYPE": "real" if self.is_real else "complex",
            "$X_VARIABLE_NAME": VariableNameFromRangeKey[max_key],
            "$X_LABEL": as_string(self.get_x_label(), raw_string=True),
        }

        for atom in range(self.n_atoms):
            replacements[f"$SPECIES_{atom}"] = as_string(self.get_species(atom))
            replacements[f"$QUANTUM_NUMBERS_{atom}"] = dict_to_repl(self.get_quantum_numbers(atom))
            replacements[f"$QUANTUM_NUMBERS_RESTRICTIONS_{atom}"] = dict_to_repl(
                self.get_quantum_number_restrictions(atom)
            )

        replacements["$STEPS"] = str(self.steps)
        for key, values in self.ranges.items():
            replacements[f"${key.upper()}_MIN"] = str(values[0])
            replacements[f"${key.upper()}_MAX"] = str(values[-1])
            if values[0] == values[-1]:
                replacements[f"${key.upper()}_VALUE"] = str(values[0])

        replacements["$DIAGONALIZE_KWARGS"] = dict_to_repl(self.diagonalize_kwargs)

        if self.diagonalize_relative_energy_range is not None:
            r_energy = self.diagonalize_relative_energy_range
            replacements["$DIAGONALIZE_ENERGY_RANGE_KWARGS"] = (
                f', energy_range=(ket_energy + {r_energy[0]}, ket_energy - {-r_energy[1]}), energy_unit="GHz"'
            )
        else:
            replacements["$DIAGONALIZE_ENERGY_RANGE_KWARGS"] = ""

        return replacements


@dataclass
class Results(ABC):
    energies: list["NDArray"]
    energy_offset: float
    ket_overlaps: list["NDArray"]
    state_labels: dict[int, list[str]]

    @classmethod
    def from_calculate(
        cls,
        parameters: Parameters[Any],
        system_list: Union[
            list[pi_real.SystemPair], list[pi_complex.SystemPair], list[pi_real.SystemAtom], list[pi_complex.SystemAtom]
        ],
        ket: Union[_wrapped.KetAtom, tuple[_wrapped.KetAtom, ...]],
        energy_offset: float,
    ) -> "Self":
        """Create Results object from ket, basis, and diagonalized systems."""
        energies = [system.get_eigenenergies("GHz") - energy_offset for system in system_list]
        ket_overlaps = [system.get_eigenbasis().get_overlaps(ket) for system in system_list]  # type: ignore [arg-type]

        steps_with_labels = [int(i) for i in np.linspace(0, parameters.steps - 1, parameters.number_state_labels)]
        states_dict = {i: system_list[i].get_eigenbasis().states for i in steps_with_labels}
        state_labels = {i: [s.get_label() for s in states] for i, states in states_dict.items()}

        return cls(energies, energy_offset, ket_overlaps, state_labels)


def as_string(value: str, *, raw_string: bool = False) -> str:
    string = '"' + value + '"'
    if raw_string:
        string = "r" + string
    return string


def dict_to_repl(d: dict[str, Any]) -> str:
    """Convert a dictionary to a string for replacement."""
    if not d:
        return ""
    repl = ""
    for k, v in d.items():
        if isinstance(v, str):
            repl += f", {k}={as_string(v)}"
        else:
            repl += f", {k}={v}"
    return repl
