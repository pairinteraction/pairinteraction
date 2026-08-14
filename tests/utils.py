# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

import contextlib
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Protocol

import numpy as np
import pytest

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator

    import pairinteraction as pi
    import pint
    from pairinteraction.units import NDArray


REFERENCE_PATHS = {
    "stark_map": Path(__file__).parent.parent / "data" / "reference_stark_map",
    "pair_potential": Path(__file__).parent.parent / "data" / "reference_pair_potential",
}

SHRUNK_DATABASE_PATH = (Path(__file__).parent.parent / "data" / "database").resolve()


def is_shrunk_database_used() -> bool:
    """Check whether the shrunk database that comes with the repository is used.

    The shrunk database only contains states within a narrow range of the effective principal quantum number, a
    few low-lying states, and only small values of the quantum number l_ryd. Thus, results that depend on the
    completeness of the database, like lifetimes, are far off and must not be checked against reference values.
    To check these values as well, run the tests with `pytest --database-dir "" --download-missing`.
    """
    from pairinteraction import Database

    return Path(Database.get_global_database().database_dir).resolve() == SHRUNK_DATABASE_PATH


def skip_value_check_if_shrunk_database() -> None:
    """Skip the rest of the current test if the results cannot be checked against reference values.

    Call this directly before the checks of the resulting values, see :func:`is_shrunk_database_used`. The code
    before the call is still executed, only the checks are skipped and the test is reported as skipped.
    """
    if is_shrunk_database_used():
        pytest.skip(
            "the shrunk database is used, run with `pytest --database-dir '' --download-missing` to include all checks"
        )


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


class PairinteractionModule(Protocol):
    ureg: pint.UnitRegistry
    Database: type[pi.Database]
    KetAtom: type[pi.KetAtom]
    StateAtom: type[pi.StateAtom]
    BasisAtom: type[pi.BasisAtom]
    SystemAtom: type[pi.SystemAtom]
    KetPair: type[pi.KetPair]
    BasisPair: type[pi.BasisPair]
    SystemPair: type[pi.SystemPair]
    EffectiveSystemPair: type[pi.EffectiveSystemPair]
    C3: type[pi.C3]
    C6: type[pi.C6]
    diagonalize: Callable[..., None]
