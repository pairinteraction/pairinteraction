# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from pathlib import Path
from typing import TYPE_CHECKING, Optional

import numpy as np

if TYPE_CHECKING:
    from pairinteraction.units import NDArray


REFERENCE_PATHS = {
    "stark_map": Path(__file__).parent.parent / "data" / "reference_stark_map",
    "pair_potential": Path(__file__).parent.parent / "data" / "reference_pair_potential",
}


def compare_to_reference(
    reference_path: Path,
    eigenenergies: "NDArray",
    overlaps: Optional["NDArray"] = None,
    eigenvectors: Optional["NDArray"] = None,
    kets: Optional[list[str]] = None,
) -> None:
    n_systems, n_kets = eigenenergies.shape
    np.testing.assert_allclose(eigenenergies, np.loadtxt(reference_path / "eigenenergies.txt"))

    if overlaps is not None:
        # Ensure that the overlaps sum up to one
        np.testing.assert_allclose(np.sum(overlaps, axis=1), np.ones(n_systems))
        np.testing.assert_allclose(overlaps, np.loadtxt(reference_path / "overlaps.txt"), rtol=1e-3, atol=1e-4)

    if kets is not None:
        np.testing.assert_equal(kets, np.loadtxt(reference_path / "kets.txt", dtype=str, delimiter="\t"))

    if eigenvectors is not None:
        # Because of degeneracies, checking the eigenvectors against reference data is complicated.
        # Thus, we only check their normalization and orthogonality.
        cumulative_norm = (np.array(eigenvectors) * np.array(eigenvectors).conj()).sum(axis=1)
        np.testing.assert_allclose(cumulative_norm, n_kets * np.ones(n_systems))
