# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from typing import TYPE_CHECKING, Optional, Union

import numpy as np
from attr import dataclass

from pairinteraction import (
    complex as pi_complex,
    real as pi_real,
)
from pairinteraction_gui.calculate.calculate_base import Parameters, Results
from pairinteraction_gui.worker import run_in_other_process

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction_gui.page import TwoAtomsPage

logger = logging.getLogger(__name__)


@dataclass
class ParametersTwoAtoms(Parameters["TwoAtomsPage"]):
    """Parameters for the two atoms calculation."""

    pair_delta_energy: float = np.inf
    pair_m_range: Optional[tuple[float, float]] = None
    order: int = 3

    @classmethod
    def from_page(cls, page: "TwoAtomsPage") -> "Self":
        obj = super().from_page(page)
        obj.pair_delta_energy = page.basis_config.pair_delta_energy.value(np.inf)
        obj.pair_m_range = (
            page.basis_config.pair_m_range.values() if page.basis_config.pair_m_range.isChecked() else None
        )
        obj.order = page.system_config.order.value()
        return obj

    def to_replacement_dict(self) -> dict[str, str]:
        replacements = super().to_replacement_dict()
        replacements["$MULTIPOLE_ORDER"] = str(self.order)
        replacements["$PAIR_DELTA_ENERGY"] = (
            "np.inf" if np.isinf(self.pair_delta_energy) else str(self.pair_delta_energy)
        )
        replacements["$PAIR_M_RANGE"] = str(self.pair_m_range)
        return replacements


@dataclass
class ResultsTwoAtoms(Results):
    basis_0_label: Optional[str] = None


@run_in_other_process
def calculate_two_atoms(parameters: ParametersTwoAtoms) -> ResultsTwoAtoms:
    """Calculate the energy plot for two atoms.

    This means, given a Parameters object, do the pairinteraction calculations and return an ResultsTwoAtoms object.
    """
    return _calculate_two_atoms(parameters)


def _calculate_two_atoms(parameters: ParametersTwoAtoms) -> ResultsTwoAtoms:
    """Make the unwrapped function available for testing."""
    pi = pi_real if parameters.is_real else pi_complex
    n_atoms = 2

    kets = tuple(pi.KetAtom(parameters.get_species(i), **parameters.get_quantum_numbers(i)) for i in range(n_atoms))
    bases = tuple(
        pi.BasisAtom(parameters.get_species(i), **parameters.get_quantum_number_restrictions(i)) for i in range(n_atoms)
    )

    fields = {k: v for k, v in parameters.ranges.items() if k in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]}

    basis_pair_list: Union[list[pi_real.BasisPair], list[pi_complex.BasisPair]]
    if all(v[0] == v[-1] for v in fields.values()):
        # If all fields are constant, we can only have to diagonalize one SystemAtom per atom
        # and can construct one BasisPair, which we can use for all steps
        systems = tuple(
            pi.SystemAtom(bases[i])
            .set_electric_field(parameters.get_efield(0), unit="V/cm")
            .set_magnetic_field(parameters.get_bfield(0), unit="G")
            for i in range(n_atoms)
        )
        logger.debug("Diagonalizing SystemAtoms...")
        pi.diagonalize(systems, **parameters.diagonalize_kwargs)
        logger.debug("Done diagonalizing SystemAtoms.")
        ket_pair_energy_0 = sum(systems[i].get_corresponding_energy(kets[i], "GHz") for i in range(n_atoms))
        delta_energy = parameters.pair_delta_energy
        basis_pair = pi.BasisPair(
            systems,
            energy=(ket_pair_energy_0 - delta_energy, ket_pair_energy_0 + delta_energy),
            energy_unit="GHz",
            m=parameters.pair_m_range,
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
        logger.debug("Diagonalizing SystemAtoms...")
        pi.diagonalize(systems_flattened, **parameters.diagonalize_kwargs)
        logger.debug("Done diagonalizing SystemAtoms.")
        delta_energy = parameters.pair_delta_energy
        basis_pair_list = []
        for step in range(parameters.steps):
            ket_pair_energy = sum(
                systems_list[step][i].get_corresponding_energy(kets[i], "GHz") for i in range(n_atoms)
            )
            basis_pair = pi.BasisPair(
                systems_list[step],
                energy=(ket_pair_energy - delta_energy, ket_pair_energy + delta_energy),
                energy_unit="GHz",
                m=parameters.pair_m_range,
            )
            basis_pair_list.append(basis_pair)
        ket_pair_energy_0 = sum(systems_list[-1][i].get_corresponding_energy(kets[i], "GHz") for i in range(n_atoms))

    system_pair_list: Union[list[pi_real.SystemPair], list[pi_complex.SystemPair]] = []
    for step in range(parameters.steps):
        system = pi.SystemPair(basis_pair_list[step])
        system.set_interaction_order(parameters.order)
        if "Distance" in parameters.ranges:
            distance = parameters.ranges["Distance"][step]
            angle: float = 0
            if "Angle" in parameters.ranges:
                angle = parameters.ranges["Angle"][step]
            system.set_distance(distance, angle, unit="micrometer")
        system_pair_list.append(system)

    logger.debug("Diagonalizing SystemPairs...")
    pi.diagonalize(
        system_pair_list,
        **parameters.diagonalize_kwargs,
        **parameters.get_diagonalize_energy_range(ket_pair_energy_0),
    )
    logger.debug("Done diagonalizing SystemPairs.")

    results = ResultsTwoAtoms.from_calculate(parameters, system_pair_list, kets, ket_pair_energy_0)
    results.basis_0_label = (
        str(basis_pair_list[-1]) + f"\n  ⇒ Basis consists of {basis_pair_list[-1].number_of_kets} kets"
    )

    return results
