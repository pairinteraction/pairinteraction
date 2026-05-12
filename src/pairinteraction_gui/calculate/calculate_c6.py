# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

import pairinteraction as pi
from pairinteraction.units import QuantityArray
from pairinteraction_gui.calculate.calculate_two_atoms import ParametersTwoAtoms

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from pairinteraction.ket import KetPair


logger = logging.getLogger(__name__)


class ParametersC6(ParametersTwoAtoms):
    """Parameters for the C6 calculation."""


@dataclass
class ResultsC6:
    """Results of the C6 calculation.

    All energies are in GHz, contributions to the C6 coefficient are in GHz·μm⁶
    so that ``contributions_c6.sum()`` reproduces ``c6``.
    """

    c6_obj: pi.C6
    c6: float
    kets: list[KetPair]
    gaps_ghz: NDArray[np.float64]
    couplings_ghz: NDArray[np.float64]
    contributions_c6: NDArray[np.float64]


def calculate_c6(parameters: ParametersC6) -> ResultsC6:
    """Calculate the C6 coefficient and per-intermediate-state contributions."""
    ket1 = pi.KetAtom(parameters.species[0], **parameters.quantum_numbers[0])
    ket2 = pi.KetAtom(parameters.species[1], **parameters.quantum_numbers[1])
    c6_obj = pi.C6(ket1, ket2)
    c6_obj.set_angle(parameters.get_angle(0), unit="degree")
    c6_obj.set_interaction_order(parameters.order)

    basis1 = pi.BasisAtom(parameters.species[0], **parameters.quantum_number_restrictions[0])
    basis2 = pi.BasisAtom(parameters.species[1], **parameters.quantum_number_restrictions[1])
    basis_atoms = (basis1, basis2)

    system_atoms = tuple(
        pi.SystemAtom(basis)
        .set_electric_field(parameters.get_efield(0), unit="V/cm")
        .set_magnetic_field(parameters.get_bfield(0), unit="G")
        .set_diamagnetism_enabled(parameters.diamagnetism_enabled)
        for basis in basis_atoms
    )

    c6_obj.system_atoms = system_atoms

    if np.isfinite(parameters.pair_delta_energy):
        c6_obj.create_basis_pair(delta_energy=parameters.pair_delta_energy, delta_energy_unit="GHz")

    c6 = c6_obj.get(unit="planck_constant * GHz * micrometer^6")

    ket_pairs, gaps_pint, couplings_pint = c6_obj.get_kets_with_contributions()
    gaps_ghz = QuantityArray.from_pint(gaps_pint, "energy").to_unit("GHz")
    couplings_ghz = QuantityArray.from_pint(couplings_pint, "energy").to_unit("GHz")

    # C_6 = sum_k |coupling_k|^2 / gap_k * r^6 (see C6.get in src/pairinteraction/perturbative/c6.py).
    # Express each intermediate state's contribution in GHz·μm⁶ so the bars
    # sum to the displayed C6 value.
    distance_um = c6_obj.system_pair.get_distance("micrometer")
    with np.errstate(divide="ignore", invalid="ignore"):
        contributions_c6 = np.where(gaps_ghz != 0, np.abs(couplings_ghz) ** 2 / gaps_ghz * distance_um**6, np.inf)

    return ResultsC6(
        c6_obj=c6_obj,
        c6=c6,
        kets=ket_pairs,
        gaps_ghz=gaps_ghz,
        couplings_ghz=couplings_ghz,
        contributions_c6=contributions_c6,
    )
