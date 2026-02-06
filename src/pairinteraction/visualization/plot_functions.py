# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from pairinteraction.ket.ket_pair import KetPair, is_ket_atom_tuple
from pairinteraction.visualization.colormaps import alphamagma

if TYPE_CHECKING:
    from collections.abc import Sequence

    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    from pairinteraction.ket.ket_pair import KetPairLike
    from pairinteraction.system.system_pair import SystemPair
    from pairinteraction.units import NDArray


def plot_pair_potential(
    system_pairs: Sequence[SystemPair],
    ket_pair_of_interest: KetPairLike,
    ax: Axes | None = None,
) -> tuple[Figure, Axes]:
    """Plot the pair potential of a sequence of SystemPair objects.

    Args:
        system_pairs: A sequence of SystemPair objects to plot.
        ket_pair_of_interest: The KetPair (or tuple of KetAtom) for whicht the overlaps will be plotted.
        ax: An optional matplotlib Axes object to plot on.
            If None, a new figure and axes will be created.

    """
    distances: NDArray = np.array([system.get_distance(unit="micrometer") for system in system_pairs])

    if is_ket_atom_tuple(ket_pair_of_interest):
        energy_of_interest = sum(ket.get_energy(unit="GHz") for ket in ket_pair_of_interest)
    elif isinstance(ket_pair_of_interest, KetPair):
        energy_of_interest = ket_pair_of_interest.get_energy(unit="GHz")
    else:
        raise TypeError("ket_pair_of_interest must be a KetPair or a tuple of KetAtom.")

    eigenenergies = [system.get_eigenenergies(unit="GHz") - energy_of_interest for system in system_pairs]
    overlaps = [system.get_eigenbasis().get_overlaps(ket_pair_of_interest) for system in system_pairs]

    if ax is None:
        owns_ax = True
        fig, ax = plt.subplots()
    else:
        owns_ax = False
        fig = ax.figure  # type: ignore[assignment]

    lw = mpl.rcParams["lines.linewidth"]
    ax.plot(distances, np.array(eigenenergies), c="k", lw=lw / 6, zorder=-10)

    x_repeated = np.hstack([val * np.ones_like(es) for val, es in zip(distances, eigenenergies)])
    energies_flattend = np.hstack(eigenenergies)
    overlaps_flattend = np.hstack(overlaps)
    sorter = np.argsort(overlaps_flattend)

    ms = mpl.rcParams["lines.markersize"]
    scat = ax.scatter(
        x_repeated[sorter],
        energies_flattend[sorter],
        c=overlaps_flattend[sorter],
        s=ms * 2.5,
        vmin=0,
        vmax=1,
        cmap=alphamagma,
    )

    if owns_ax:
        ax.set_xlabel(r"Distance ($\mu$m)")
        ax.set_ylabel(r"Energy (h$\cdot$GHz)")
        fig.colorbar(scat, ax=ax, label="Overlap with state of interest")
        fig.tight_layout()

    return fig, ax
