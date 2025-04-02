from typing import TYPE_CHECKING, Any, Union

from attr import dataclass

from pairinteraction import (
    _wrapped,
    complex as pi_complex,
    real as pi_real,
)

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix
    from typing_extensions import Self

    from pairinteraction.units import NDArray
    from pairinteraction_gui.page import OneAtomPage

# FIXME: having all kwargs dictionaries being Any is a hacky solution, it would be nice to use TypedDict in the future


@dataclass
class ParametersOneAtom:
    species: str
    quantum_numbers: dict[str, float]
    quantum_number_deltas: dict[str, float]
    ranges: dict[str, list[float]]
    diagonalize_kwargs: dict[str, str]
    diagonalize_relative_energy_range: Union[tuple[float, float], None]

    @classmethod
    def from_one_atom_page(cls, page: "OneAtomPage") -> "Self":
        """Create ParametersOneAtom object from one_atom_page."""
        species = page.ket_config.get_species()
        quantum_numbers = page.ket_config.get_quantum_numbers()

        quantum_number_deltas = page.basis_config.get_quantum_number_deltas()

        ranges = page.system_atom_config.get_ranges_dict()

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
            ranges,
            diagonalize_kwargs,
            diagonalize_relative_energy_range,
        )

    def is_real(self) -> bool:
        """Check if the parameters are real."""
        return all(e == 0 for e in self.ranges.get("Ey", [0])) and all(b == 0 for b in self.ranges.get("By", [0]))

    def get_efield(self, step: int) -> list[float]:
        """Return the electric field for the given step."""
        return [self.ranges[key][step] if key in self.ranges else 0 for key in ["Ex", "Ey", "Ez"]]

    def get_bfield(self, step: int) -> list[float]:
        """Return the magnetic field for the given step."""
        return [self.ranges[key][step] if key in self.ranges else 0 for key in ["Bx", "By", "Bz"]]

    @property
    def steps(self) -> int:
        """Return the number of steps."""
        steps = len(next(iter(self.ranges.values())))
        assert all(len(v) == steps for v in self.ranges.values()), "All ranges must have the same number of steps"
        return steps

    def get_quantum_numbers(self) -> dict[str, Any]:
        """Return the quantum numbers for the given ket."""
        return self.quantum_numbers

    def get_quantum_number_restrictions(self, ket: _wrapped.KetAtom) -> dict[str, Any]:
        """Return the quantum number restrictions for the given ket."""
        qn_restrictions: dict[str, tuple[float, float]] = {}
        for key, delta in self.quantum_number_deltas.items():
            if key in self.quantum_numbers:
                qn_restrictions[key] = (self.quantum_numbers[key] - delta, self.quantum_numbers[key] + delta)
            else:
                raise ValueError(f"Quantum number delta {key} not found in quantum numbers.")
        return qn_restrictions

    def get_diagonalize_kwargs(self, ket: _wrapped.KetAtom) -> dict[str, Any]:
        """Return the kwargs for the diagonalization function."""
        kwargs: dict[str, Any] = self.diagonalize_kwargs
        if self.diagonalize_relative_energy_range is not None:
            kwargs["energy_range"] = (
                ket.get_energy("GHz") + self.diagonalize_relative_energy_range[0],
                ket.get_energy("GHz") + self.diagonalize_relative_energy_range[1],
            )
            kwargs["energy_unit"] = "GHz"

        return kwargs

    def _get_max_diff_key(self) -> str:
        """Return the key with the maximum difference in the ranges."""
        range_diffs = {key: abs(r[-1] - r[0]) for key, r in self.ranges.items()}
        return max(range_diffs, key=range_diffs.get)

    def get_x_values(self) -> list[float]:
        """Return the x values for the plot."""
        max_key = self._get_max_diff_key()
        return self.ranges[max_key]

    def get_x_label(self) -> str:
        """Return the x values for the plot."""
        max_key = self._get_max_diff_key()

        x_label = max_key
        x_units = {"E": "V/cm", "B": "Gauss", "distance": r"$\mu$m", "angle": r"$^\circ$"}
        x_label += f" [{next((unit for key, unit in x_units.items() if max_key.startswith(key)), '')}]"

        non_constant_keys = [key for key, values in self.ranges.items() if key != max_key and values[0] != values[-1]]
        if non_constant_keys:
            x_label += f"  ({', '.join(non_constant_keys)} did also change)"

        return x_label


@dataclass
class ResultsOneAtom:
    ket_energy: float
    energies: list["NDArray"]
    ket_overlaps: list["NDArray"]
    all_overlaps: list["csr_matrix"]
    corresponding_kets_0: list[str]

    @classmethod
    def from_ket_basis_systems(
        cls,
        ket: _wrapped.KetAtom,
        basis: Union[pi_real.BasisAtom, pi_complex.BasisAtom],
        systems: Union[list[pi_real.SystemAtom], list[pi_complex.SystemAtom]],
    ) -> "Self":
        """Create ResultsOneAtom object from ket, basis, and diagonalized systems."""
        ket_energy = ket.get_energy("GHz")
        energies = [system.get_eigenenergies("GHz") - ket_energy for system in systems]
        ket_overlaps = [system.get_eigenbasis().get_overlaps(ket) for system in systems]
        all_overlaps = [system.get_eigenbasis().get_overlaps(basis) for system in systems]
        basis_0 = systems[0].get_eigenbasis()
        corresponding_kets_0 = [
            basis_0.get_corresponding_ket(i).get_label("ket") for i in range(basis_0.number_of_states)
        ]

        return cls(ket_energy, energies, ket_overlaps, all_overlaps, corresponding_kets_0)


def calculate_one_atom(parameters: ParametersOneAtom) -> ResultsOneAtom:
    """Calculate the energy plot for one atom.

    This means, given a Paramaters object, do the pairinteraction calculations and return an ResultsOneAtom object.
    """
    pi = pi_real if parameters.is_real() else pi_complex

    ket = pi.KetAtom(parameters.species, **parameters.get_quantum_numbers())
    basis = pi.BasisAtom(parameters.species, **parameters.get_quantum_number_restrictions(ket))

    system_list = [
        pi.SystemAtom(basis)
        .set_electric_field(parameters.get_efield(step), unit="V/cm")
        .set_magnetic_field(parameters.get_bfield(step), unit="G")
        for step in range(parameters.steps)
    ]

    pi.diagonalize(system_list, **parameters.get_diagonalize_kwargs(ket))

    return ResultsOneAtom.from_ket_basis_systems(ket, basis, system_list)
