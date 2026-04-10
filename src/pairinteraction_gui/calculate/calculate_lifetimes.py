# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import pairinteraction as pi
from pairinteraction_gui.calculate.calculate_base import as_string, dict_to_repl

if TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray

    from pairinteraction_gui.config.ket_config import QuantumNumbers
    from pairinteraction_gui.page import LifetimesPage

logger = logging.getLogger(__name__)


@dataclass
class ParametersLifetimes:
    """Parameters for the lifetimes calculation."""

    species: str
    quantum_numbers: QuantumNumbers
    temperature: float

    @classmethod
    def from_page(cls, page: LifetimesPage) -> ParametersLifetimes:
        """Create ParametersLifetimes from the GUI page state."""
        return cls(
            species=page.ket_config.get_species(),
            quantum_numbers=page.ket_config.get_quantum_numbers(),
            temperature=page.ket_config.get_temperature(),
        )

    def to_replacement_dict(self) -> dict[str, str]:
        """Return a dictionary with the parameters for template replacement."""
        return {
            "$SPECIES_0": as_string(self.species),
            "$QUANTUM_NUMBERS_0": dict_to_repl(self.quantum_numbers),
            "$TEMPERATURE": str(self.temperature),
        }


@dataclass
class KetData:
    """Picklable representation of a KetAtom for inter-process transfer."""

    n: int
    label: str

    def __str__(self) -> str:
        return self.label


@dataclass
class ResultsLifetimes:
    """Results for the lifetimes calculation."""

    kets_sp: list[KetData]
    transition_rates_sp: NDArray[np.float64]
    kets_bbr: list[KetData]
    transition_rates_bbr: NDArray[np.float64]
    lifetime: float


def calculate_lifetimes(parameters: ParametersLifetimes) -> ResultsLifetimes:
    """Calculate the transition rates for a given ket."""
    return _calculate_lifetimes(parameters)


def _calculate_lifetimes(parameters: ParametersLifetimes) -> ResultsLifetimes:
    """Make the unwrapped function available for testing."""
    ket = pi.KetAtom(parameters.species, **parameters.quantum_numbers)
    kets_sp, rates_sp = ket.get_spontaneous_transition_rates(unit="1/ms")
    kets_bbr, rates_bbr = ket.get_black_body_transition_rates(parameters.temperature, "K", unit="1/ms")
    lifetime = ket.get_lifetime(parameters.temperature, temperature_unit="K", unit="mus")
    return ResultsLifetimes(
        kets_sp=[KetData(n=int(k.n), label=str(k)) for k in kets_sp],
        transition_rates_sp=rates_sp,
        kets_bbr=[KetData(n=int(k.n), label=str(k)) for k in kets_bbr],
        transition_rates_bbr=rates_bbr,
        lifetime=lifetime,
    )
