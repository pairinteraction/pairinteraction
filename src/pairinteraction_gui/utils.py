# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import re
from typing import Literal

from pairinteraction.cli import download_databases
from pairinteraction_gui.worker import run_in_other_process

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
        table = next(w for w in str(err).split(" ") if "states" in w)
        self.species = table.replace("_states", "")


class NoStateFoundError(Exception):
    def __init__(self, err: ValueError) -> None:
        super().__init__(str(err))


def get_custom_error(err: Exception) -> Exception:
    """Get a custom error message based on the type of error."""
    if isinstance(err, RuntimeError) and "Table" in str(err) and "not found" in str(err):
        return DatabaseMissingError(err)
    if isinstance(err, ValueError) and ("No state found" in str(err) or "quantum number m must be" in str(err)):
        return NoStateFoundError(err)
    return err


@run_in_other_process
def download_databases_mp(species: list[str]) -> None:
    """Download the databases in a separate process."""
    download_databases(species)


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
