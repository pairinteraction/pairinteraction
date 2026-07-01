# SPDX-FileCopyrightText: 2026 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

QUANTUM_NUMBER_L_LABELS = ("S", "P", "D", "F", "G", "H")


def format_half_integer(value: float) -> str:
    """Format an integer or half-integer quantum number as e.g. "1" or "1/2"."""
    if value == round(value):
        return f"{value:.0f}"
    if 2 * value == round(2 * value):
        return f"{2 * value:.0f}/2"
    logger.warning("Value %.5f is not an integer or half-integer.", value)
    return f"{value:.1f}"


def get_l_label(l: float) -> str:
    """Get the label for an integer quantum number l."""
    if l == round(l):
        if l < len(QUANTUM_NUMBER_L_LABELS):
            return QUANTUM_NUMBER_L_LABELS[round(l)]
        return f"{round(l):d}"
    logger.warning("Value of l %.5f is not an integer.", l)
    return f"{l:.1f}"
