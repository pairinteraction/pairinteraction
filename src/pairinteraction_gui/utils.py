# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import re
from typing import TYPE_CHECKING, Literal

import numpy as np

if TYPE_CHECKING:
    from pairinteraction.units import NDArray

AVAILABLE_SPECIES = [
    "Rb",
    "Li",
    "Na",
    "K",
    "Cs",
    "Sr88_singlet",
    "Sr88_triplet",
    "Sr87_mqdt",
    "Sr88_mqdt",
    "Yb171_mqdt",
    "Yb173_mqdt",
    "Yb174_mqdt",
]
SpeciesTypes = Literal["sqdt_duplet", "sqdt_singlet", "sqdt_triplet", "mqdt_halfint", "mqdt_int"]


class DatabaseMissingError(Exception):
    def __init__(self, err: RuntimeError) -> None:
        super().__init__(str(err))


class NoStateFoundError(Exception):
    def __init__(self, err: ValueError) -> None:
        super().__init__(str(err))


def get_custom_error(err: Exception) -> Exception:
    """Get a custom error message based on the type of error."""
    if isinstance(err, RuntimeError) and "No tables found for" in str(err):
        return DatabaseMissingError(err)
    if isinstance(err, ValueError) and ("No state found" in str(err) or "quantum number m must be" in str(err)):
        return NoStateFoundError(err)
    return err


def get_species_type(species: str) -> SpeciesTypes:
    """Return the species type based on the species name of the ... atom."""
    if "mqdt" in species:
        match = re.search(r"\d+", species)
        if match:
            if int(match.group()) % 2 == 0:
                return "mqdt_int"
            return "mqdt_halfint"
        raise ValueError(f"Invalid species name: {species}")
    if "singlet" in species:
        return "sqdt_singlet"
    if "triplet" in species:
        return "sqdt_triplet"
    return "sqdt_duplet"


def label_to_object_name(label: str) -> str:
    """Convert a display label to an human readable object name (and QSettings key)."""
    label = label.lower().strip()
    label = re.sub(r"\s+", "_", label)
    label = re.sub(r"[\\/]+", "_", label)
    return label.replace("δ", "delta")


def distance_angle_to_xz_vector(distance: float, angle_degree: float) -> NDArray:
    """Convert a distance-angle description to an xz-plane vector.

    0° corresponds to the positive z-axis and 90° to the positive x-axis.
    """
    angle_rad = np.deg2rad(angle_degree)
    return np.array([np.sin(angle_rad) * distance, 0.0, np.cos(angle_rad) * distance])
