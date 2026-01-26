# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

from pairinteraction.ket.ket_pair import KetPair, is_ket_atom_tuple
from pairinteraction.visualization.colormaps import alphamagma

if TYPE_CHECKING:
    from collections.abc import Sequence

    from pairinteraction.ket.ket_pair import KetPairLike
    from pairinteraction.system.system_pair import SystemPair
    from pairinteraction.units import NDArray


def plot_pair_potential(
    system_pairs: Sequence[SystemPair], ket_pair_of_interest: KetPairLike
) -> tuple[plt.Figure, plt.Axes]:
    """Plot the pair potential of a sequence of SystemPair objects."""
    distances: NDArray = np.array([system.get_distance(unit="micrometer") for system in system_pairs])

    if is_ket_atom_tuple(ket_pair_of_interest):
        energy_of_interest = sum(ket.get_energy(unit="GHz") for ket in ket_pair_of_interest)
    elif isinstance(ket_pair_of_interest, KetPair):
        energy_of_interest = ket_pair_of_interest.get_energy(unit="GHz")
    else:
        raise TypeError("ket_pair_of_interest must be a KetPair or a tuple of KetAtom.")

    eigenenergies = [system.get_eigenenergies(unit="GHz") - energy_of_interest for system in system_pairs]
    overlaps = [system.get_eigenbasis().get_overlaps(ket_pair_of_interest) for system in system_pairs]

    fig, ax = plt.subplots()

    ax.set_xlabel(r"Distance [$\mu$m]")
    ax.set_ylabel("Energy [GHz]")

    ax.plot(distances, np.array(eigenenergies), c="k", lw=0.25, zorder=-10)

    x_repeated = np.hstack([val * np.ones_like(es) for val, es in zip(distances, eigenenergies)])
    energies_flattend = np.hstack(eigenenergies)
    overlaps_flattend = np.hstack(overlaps)
    sorter = np.argsort(overlaps_flattend)

    scat = ax.scatter(
        x_repeated[sorter],
        energies_flattend[sorter],
        c=overlaps_flattend[sorter],
        s=15,
        vmin=0,
        vmax=1,
        cmap=alphamagma,
    )

    fig.colorbar(scat, ax=ax, label="Overlap with state of interest")

    return fig, ax
