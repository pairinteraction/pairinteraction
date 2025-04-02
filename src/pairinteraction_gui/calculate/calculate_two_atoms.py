from typing import TYPE_CHECKING, Any, Optional, Union

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
    from pairinteraction_gui.page import TwoAtomsPage

# FIXME: having all kwargs dictionaries being Any is a hacky solution, it would be nice to use TypedDict in the future


@dataclass
class ParametersTwoAtoms:
    species: tuple[str, ...]
    quantum_numbers: tuple[dict[str, float], ...]
    quantum_number_deltas: tuple[dict[str, float], ...]
    pair_basis_energy_delta: float
    ranges: dict[RangesKeys, list[float]]
    order: int
    diagonalize_kwargs: dict[str, str]
    diagonalize_relative_energy_range: Union[tuple[float, float], None]

    @classmethod
    def from_two_atoms_page(cls, page: "TwoAtomsPage") -> "Self":
        """Create ParametersTwoAtoms object from one_atom_page."""
        n_atoms = 2

        species = tuple(page.ket_config.get_species(atom) for atom in range(n_atoms))
        quantum_numbers = tuple(page.ket_config.get_quantum_numbers(atom) for atom in range(n_atoms))

        quantum_number_deltas = tuple(page.basis_config.get_quantum_number_deltas(atom) for atom in range(n_atoms))
        pair_basis_energy_delta = page.basis_config.delta_pair_energy.value()

        ranges = page.system_pair_config.get_ranges_dict()

        order = page.system_pair_config.order.value()

        diagonalize_kwargs = {}
        if page.plotwidget.fast_mode.isChecked():
            diagonalize_kwargs["diagonalizer"] = "lapacke_evr"
            diagonalize_kwargs["float_type"] = "float32"

        diagonalize_relative_energy_range = None
        if page.plotwidget.energy_range.isChecked():
            diagonalize_relative_energy_range = page.plotwidget.energy_range.values()

        return cls(
            species,
            quantum_numbers,
            quantum_number_deltas,
            pair_basis_energy_delta,
            ranges,
            order,
            diagonalize_kwargs,
            diagonalize_relative_energy_range,
        )

    def is_real(self) -> bool:
        """Check if the parameters are real."""
        return all(e == 0 for e in self.ranges.get("Ey", [0])) and all(b == 0 for b in self.ranges.get("By", [0]))

    def get_efield(self, step: int) -> list[float]:
        """Return the electric field for the given step."""
        efield_keys: list[RangesKeys] = ["Ex", "Ey", "Ez"]
        return [self.ranges[key][step] if key in self.ranges else 0 for key in efield_keys]

    def get_bfield(self, step: int) -> list[float]:
        """Return the magnetic field for the given step."""
        bfield_keys: list[RangesKeys] = ["Bx", "By", "Bz"]
        return [self.ranges[key][step] if key in self.ranges else 0 for key in bfield_keys]

    @property
    def steps(self) -> int:
        """Return the number of steps."""
        steps = len(next(iter(self.ranges.values())))
        assert all(len(v) == steps for v in self.ranges.values()), "All ranges must have the same number of steps"
        return steps

    def get_species(self, atom: Optional[int] = None) -> str:
        """Return the species for the given ket."""
        if atom is None:
            assert len(self.species) == 1, "Please specify the atom index when using multiple atoms."
            atom = 0
        return self.species[atom]

    def get_quantum_numbers(self, atom: Optional[int] = None) -> dict[str, Any]:
        """Return the quantum numbers for the given ket."""
        if atom is None:
            assert len(self.quantum_numbers) == 1, "Please specify the atom index when using multiple atoms."
            atom = 0
        return self.quantum_numbers[atom]

    def get_quantum_number_restrictions(self, ket: _wrapped.KetAtom, atom: Optional[int] = None) -> dict[str, Any]:
        """Return the quantum number restrictions for the given ket."""
        if atom is None:
            assert len(self.quantum_number_deltas) == 1, "Please specify the atom index when using multiple atoms."
            atom = 0

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

    def get_diagonalize_kwargs(self) -> dict[str, Any]:
        """Return the kwargs for the diagonalization function."""
        return self.diagonalize_kwargs

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

    def _get_ranges_max_diff_key(self) -> RangesKeys:
        """Return the key with the maximum difference in the ranges."""
        range_diffs: dict[RangesKeys, float] = {key: abs(r[-1] - r[0]) for key, r in self.ranges.items()}
        return max(range_diffs, key=lambda x: range_diffs.get(x, -1))

    def get_x_values(self) -> list[float]:
        """Return the x values for the plot."""
        max_key = self._get_ranges_max_diff_key()
        return self.ranges[max_key]

    def get_x_label(self) -> str:
        """Return the x values for the plot."""
        max_key = self._get_ranges_max_diff_key()

        x_label: str = max_key
        units = {
            "Ex": "V/cm",
            "Ey": "V/cm",
            "Ez": "V/cm",
            "Bx": "Gauss",
            "By": "Gauss",
            "Bz": "Gauss",
            "Distance": r"$\mu$m",
            "Angle": r"$^\circ$",
        }
        x_label += f" [{units[max_key]}]"

        non_constant_keys = [key for key, values in self.ranges.items() if key != max_key and values[0] != values[-1]]
        if non_constant_keys:
            x_label += f"  ({', '.join(non_constant_keys)} did also change)"

        return x_label


@dataclass
class ResultsTwoAtoms:
    ket_pair_energy_0: float
    energies: list["NDArray"]
    ket_overlaps: list["NDArray"]
    corresponding_kets_0: list[str]
    basis_0_label: str

    @classmethod
    def from_calculate(
        cls,
        kets: tuple[_wrapped.KetAtom, ...],
        basis_pair_list: Union[list[pi_real.BasisPair], list[pi_complex.BasisPair]],
        system_pair_list: Union[list[pi_real.SystemPair], list[pi_complex.SystemPair]],
        ket_pair_energy_0: float,
    ) -> "Self":
        """Create ResultsTwoAtoms object from ket, basis, and diagonalized systems."""
        energies = [system.get_eigenenergies("GHz") - ket_pair_energy_0 for system in system_pair_list]
        ket_overlaps = [system.get_eigenbasis().get_overlaps(kets) for system in system_pair_list]
        basis_0 = system_pair_list[-1].get_eigenbasis()
        corresponding_kets_0 = [
            basis_0.get_corresponding_ket(i).get_label("ket") for i in range(basis_0.number_of_states)
        ]

        basis_0_label = str(basis_pair_list[-1]) + f"\n  ⇒ Basis consists of {basis_pair_list[-1].number_of_kets} kets"

        return cls(ket_pair_energy_0, energies, ket_overlaps, corresponding_kets_0, basis_0_label)


def calculate_two_atoms(parameters: ParametersTwoAtoms) -> ResultsTwoAtoms:
    """Calculate the energy plot for one atom.

    This means, given a Paramaters object, do the pairinteraction calculations and return an ResultsTwoAtoms object.
    """
    pi = pi_real if parameters.is_real() else pi_complex
    n_atoms = 2

    kets = tuple(pi.KetAtom(parameters.get_species(i), **parameters.get_quantum_numbers(i)) for i in range(n_atoms))
    bases = tuple(
        pi.BasisAtom(parameters.get_species(i), **parameters.get_quantum_number_restrictions(kets[i], i))
        for i in range(n_atoms)
    )

    fields = {k: v for k, v in parameters.ranges.items() if k in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]}

    # systems_list: Union[list[tuple[pi_real.SystemAtom, ...]], list[tuple[pi_complex.SystemAtom, ...]]]
    systems_list: list[tuple[Union[pi_real.SystemAtom, pi_complex.SystemAtom], ...]]
    basis_pair_list: Union[list[pi_real.BasisPair], list[pi_complex.BasisPair]]
    # basis_pair_list: list[Union[pi_real.BasisPair, pi_complex.BasisPair]]
    if all(v[0] == v[-1] for v in fields.values()):
        # If all fields are constant, we can only have to diagonalize one SystemAtom per atom
        # and can construct one BasisPair, which we can use for all steps
        systems = tuple(
            pi.SystemAtom(bases[i])
            .set_electric_field(parameters.get_efield(0), unit="V/cm")
            .set_magnetic_field(parameters.get_bfield(0), unit="G")
            for i in range(n_atoms)
        )
        pi.diagonalize(systems, **parameters.get_diagonalize_kwargs())
        ket_pair_energy_0 = sum(systems[i].get_corresponding_energy(kets[i], "GHz") for i in range(n_atoms))
        delta_energy = parameters.pair_basis_energy_delta
        basis_pair = pi.BasisPair(
            systems,
            energy=(ket_pair_energy_0 - delta_energy, ket_pair_energy_0 + delta_energy),
            energy_unit="GHz",
        )
        # not very elegant, but works (note that importantly this does not copy the basis_pair objects)
        basis_pair_list = parameters.steps * [basis_pair]
    else:
        # Otherwise, we have to diagonalize one SystemAtom per atom and per step
        # and construct one BasisPair per step
        systems_list = []
        for step in range(parameters.steps):
            systems = tuple(
                pi.SystemAtom(bases[i])
                .set_electric_field(parameters.get_efield(step), unit="V/cm")
                .set_magnetic_field(parameters.get_bfield(step), unit="G")
                for i in range(n_atoms)
            )
            systems_list.append(systems)
        systems_flattened = [system for systems in systems_list for system in systems]
        pi.diagonalize(systems_flattened, **parameters.get_diagonalize_kwargs())
        delta_energy = parameters.pair_basis_energy_delta
        basis_pair_list = []
        for step in range(parameters.steps):
            ket_pair_energy = sum(
                systems_list[step][i].get_corresponding_energy(kets[i], "GHz") for i in range(n_atoms)
            )
            basis_pair = pi.BasisPair(
                systems_list[step],
                energy=(ket_pair_energy - delta_energy, ket_pair_energy + delta_energy),
                energy_unit="GHz",
            )
            basis_pair_list.append(basis_pair)
        ket_pair_energy_0 = sum(systems_list[-1][i].get_corresponding_energy(kets[i], "GHz") for i in range(n_atoms))

    system_pair_list: Union[list[pi_real.SystemPair], list[pi_complex.SystemPair]] = []
    for step in range(parameters.steps):
        system = pi.SystemPair(basis_pair_list[step])
        system.set_order(parameters.order)
        if "Distance" in parameters.ranges:
            distance = parameters.ranges["Distance"][step]
            angle: float = 0
            if "Angle" in parameters.ranges:
                angle = parameters.ranges["Angle"][step]
            system.set_distance(distance, angle, unit="micrometer")
        system_pair_list.append(system)

    pi.diagonalize(
        system_pair_list,
        **parameters.get_diagonalize_kwargs(),
        **parameters.get_diagonalize_energy_range(ket_pair_energy_0),
    )

    return ResultsTwoAtoms.from_calculate(kets, basis_pair_list, system_pair_list, ket_pair_energy_0)
