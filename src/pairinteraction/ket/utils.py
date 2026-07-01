# SPDX-FileCopyrightText: 2026 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

QUANTUM_NUMBER_L_LABELS = ("S", "P", "D", "F", "G", "H")


def format_half_integer(value: float) -> str:
    """Format an integer or half-integer quantum number as e.g. "1" or "1/2"."""
    if value == round(value):
        return f"{value:.0f}"
    if 2 * value == round(2 * value):
        return f"{2 * value:.0f}/2"
    raise ValueError(f"Expected an integer or half-integer value, but got {value}.")
