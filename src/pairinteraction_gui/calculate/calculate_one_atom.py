# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING

from attr import dataclass

from pairinteraction import (
    complex as pi_complex,
    real as pi_real,
)
from pairinteraction_gui.calculate.calculate_base import Parameters, Results

if TYPE_CHECKING:
    from pairinteraction_gui.page import OneAtomPage  # noqa: F401  # related to ruff extend-generics


@dataclass
class ParametersOneAtom(Parameters["OneAtomPage"]):
    """Parameters for the one atom calculation."""


@dataclass
class ResultsOneAtom(Results):
    """Results for the one atom calculation."""


def calculate_one_atom(parameters: ParametersOneAtom) -> ResultsOneAtom:
    """Calculate the energy plot for one atom.

    This means, given a Paramaters object, do the pairinteraction calculations and return an ResultsOneAtom object.
    """
    pi = pi_real if parameters.is_real else pi_complex

    ket = pi.KetAtom(parameters.get_species(), **parameters.get_quantum_numbers())
    ket_energy = ket.get_energy("GHz")
    basis = pi.BasisAtom(parameters.get_species(), **parameters.get_quantum_number_restrictions())

    system_list = [
        pi.SystemAtom(basis)
        .set_electric_field(parameters.get_efield(step), unit="V/cm")
        .set_magnetic_field(parameters.get_bfield(step), unit="G")
        for step in range(parameters.steps)
    ]

    pi.diagonalize(
        system_list,
        **parameters.diagonalize_kwargs,
        **parameters.get_diagonalize_energy_range(ket_energy),
    )

    return ResultsOneAtom.from_calculate(system_list, ket, ket_energy)
