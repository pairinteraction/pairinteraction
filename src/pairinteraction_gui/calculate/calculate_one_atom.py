# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from typing import TYPE_CHECKING

from attr import dataclass

import pairinteraction as pi_complex
import pairinteraction.real as pi_real
from pairinteraction_gui.calculate.calculate_base import Parameters, Results

if TYPE_CHECKING:
    from pairinteraction_gui.page import OneAtomPage  # noqa: F401  # related to ruff extend-generics

logger = logging.getLogger(__name__)


@dataclass
class ParametersOneAtom(Parameters["OneAtomPage"]):
    """Parameters for the single-atom calculation."""


@dataclass
class ResultsOneAtom(Results):
    """Results for the single-atom calculation."""


def calculate_one_atom(parameters: ParametersOneAtom) -> ResultsOneAtom:
    """Calculate the energy plot for one atom.

    This means, given a Parameters object, do the PairInteraction calculations and return an ResultsOneAtom object.
    """
    return _calculate_one_atom(parameters)


def _calculate_one_atom(parameters: ParametersOneAtom) -> ResultsOneAtom:
    """Make the unwrapped function available for testing."""
    pi = pi_real if parameters.is_real else pi_complex

    ket = pi.KetAtom(parameters.get_species(), **parameters.get_quantum_numbers())
    ket_energy = ket.get_energy("GHz")
    basis = pi.BasisAtom(parameters.get_species(), **parameters.get_quantum_number_restrictions())

    system_list = [
        pi.SystemAtom(basis)
        .set_electric_field(parameters.get_efield(step), unit="V/cm")
        .set_magnetic_field(parameters.get_bfield(step), unit="G")
        .set_diamagnetism_enabled(parameters.diamagnetism_enabled)
        .set_ion_interaction_order(parameters.ion_interaction_order)
        .set_ion_charge(parameters.ion_charge, unit="e")
        .set_ion_distance_vector(parameters.get_ion_distance_vector(step), unit="micrometer")
        for step in range(parameters.steps)
    ]

    logger.debug("Diagonalizing SystemAtoms...")
    pi.diagonalize(
        system_list,
        **parameters.diagonalize_kwargs,
        **parameters.get_diagonalize_energy_range_kwargs(ket_energy),
    )
    logger.debug("Done diagonalizing SystemAtoms.")
    return ResultsOneAtom.from_calculate(parameters, system_list, ket, ket_energy)
