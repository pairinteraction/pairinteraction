# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import contextlib
import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterator

    from pairinteraction.units import NDArray


REFERENCE_PATHS = {
    "stark_map": Path(__file__).parent.parent / "data" / "reference_stark_map",
    "pair_potential": Path(__file__).parent.parent / "data" / "reference_pair_potential",
}


def compare_eigensystem_to_reference(
    reference_path: Path,
    eigenenergies: NDArray,
    overlaps: NDArray | None = None,
    eigenvectors: NDArray | None = None,
    kets: list[str] | None = None,
) -> None:
    n_systems, n_kets = eigenenergies.shape
    np.testing.assert_allclose(eigenenergies, np.loadtxt(reference_path / "eigenenergies.txt"))

    if overlaps is not None:
        # Ensure that the overlaps sum up to one
        np.testing.assert_allclose(np.sum(overlaps, axis=1), np.ones(n_systems))
        np.testing.assert_allclose(overlaps, np.loadtxt(reference_path / "overlaps.txt"), atol=1e-8)

    if kets is not None:
        np.testing.assert_equal(kets, np.loadtxt(reference_path / "kets.txt", dtype=str, delimiter="\t"))

    if eigenvectors is not None:
        # Because of degeneracies, checking the eigenvectors against reference data is complicated.
        # Thus, we only check their normalization and orthogonality.
        cumulative_norm = (np.array(eigenvectors) * np.array(eigenvectors).conj()).sum(axis=1)
        np.testing.assert_allclose(cumulative_norm, n_kets * np.ones(n_systems))


@contextlib.contextmanager
def no_log_propagation(logger: logging.Logger | str) -> Iterator[None]:
    """Context manager to temporarily disable log propagation for a given logger."""
    if isinstance(logger, str):
        logger = logging.getLogger(logger)
    old_value = logger.propagate
    try:
        logger.propagate = False
        yield
    finally:
        logger.propagate = old_value
