# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from typing import TYPE_CHECKING

import numpy as np
from attr import dataclass

import pairinteraction as pi_complex
import pairinteraction.real as pi_real
from pairinteraction_gui.calculate.calculate_base import Parameters, Results

if TYPE_CHECKING:
    from pairinteraction_gui.page import OneAtomPage  # noqa: F401  # related to ruff extend-generics

logger = logging.getLogger(__name__)


def _get_surface_position(distance_to_surface: float, surface_angle_degree: float) -> np.ndarray:
    """Return a point whose perpendicular distance to the plane matches the GUI input."""
    surface_angle_rad = np.deg2rad(surface_angle_degree)
    return distance_to_surface * np.array([np.sin(surface_angle_rad), 0.0, np.cos(surface_angle_rad)])


@dataclass
class ParametersOneAtom(Parameters["OneAtomPage"]):
    """Parameters for the single-atom calculation."""

    def to_replacement_dict(self) -> dict[str, str]:
        replacements = super().to_replacement_dict()
        if "DistanceToSurface" not in self.ranges:
            replacements["$GREEN_TENSOR_SURFACE_RANGE"] = ""
            replacements["$GREEN_TENSOR_SURFACE_SETUP"] = ""
            replacements["$GREEN_TENSOR_SURFACE_APPLY"] = ""
        else:
            surface_angle_values = self.ranges.get("SurfaceAngle", [0.0] * self.steps)
            replacements["$GREEN_TENSOR_SURFACE_RANGE"] = (
                f"distance_to_surface = np.linspace({self.ranges['DistanceToSurface'][0]}, "
                f"{self.ranges['DistanceToSurface'][-1]}, steps)\n"
                f"surface_angle = np.linspace({surface_angle_values[0]}, {surface_angle_values[-1]}, steps)"
            )
            replacements["$GREEN_TENSOR_SURFACE_SETUP"] = (
                "    surface_angle_rad = np.deg2rad(surface_angle[i])\n"
                "    surface_position = distance_to_surface[i] * np.array(\n"
                "        [np.sin(surface_angle_rad), 0, np.cos(surface_angle_rad)]\n"
                "    )\n"
                "    gt = pi.green_tensor.GreenTensorSurface(\n"
                "        surface_position,\n"
                "        surface_position,\n"
                "        point_on_plane=[0, 0, 0],\n"
                "        normal=[np.sin(surface_angle_rad), 0, np.cos(surface_angle_rad)],\n"
                '        unit="micrometer",\n'
                "        static_limit=True,\n"
                "    )\n"
                "    gt.without_vacuum_contribution = True\n"
            )
            replacements["$GREEN_TENSOR_SURFACE_APPLY"] = "    system.set_green_tensor(gt)\n"
        return replacements


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

    system_list = []
    for step in range(parameters.steps):
        system = (
            pi.SystemAtom(basis)
            .set_electric_field(parameters.get_efield(step), unit="V/cm")
            .set_magnetic_field(parameters.get_bfield(step), unit="G")
            .set_diamagnetism_enabled(parameters.diamagnetism_enabled)
        )
        distance_to_surface = parameters.get_distance_to_surface(step)
        if distance_to_surface is not None:
            surface_angle = parameters.get_surface_angle(step)
            surface_angle_rad = np.deg2rad(surface_angle)
            surface_position = _get_surface_position(distance_to_surface, surface_angle)
            gt = pi.green_tensor.GreenTensorSurface(
                surface_position,
                surface_position,
                point_on_plane=[0, 0, 0],
                normal=[np.sin(surface_angle_rad), 0, np.cos(surface_angle_rad)],
                unit="micrometer",
                static_limit=True,
            )
            gt.without_vacuum_contribution = True
            system.set_green_tensor(gt)
        system_list.append(system)

    logger.debug("Diagonalizing SystemAtoms...")
    pi.diagonalize(
        system_list,
        **parameters.diagonalize_kwargs,
        **parameters.get_diagonalize_energy_range_kwargs(ket_energy),
    )
    logger.debug("Done diagonalizing SystemAtoms.")
    return ResultsOneAtom.from_calculate(parameters, system_list, ket, ket_energy)
