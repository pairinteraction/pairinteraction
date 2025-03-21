# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

if TYPE_CHECKING:
    from matplotlib import colors


def _colormap_add_transparency(cmap: "colors.Colormap", min_value: float = 1e-4) -> "colors.ListedColormap":
    """Add transparency to a colormap.

    This function returns a copy of the given colormap to include an alpha channel
    that represents transparency based on the specified minimum overlap.

    Args:
        cmap: The colormap to modify.
        min_value: Reference value, at which the alpha channel should be 0.
            This also shifts the alpha channel for larger values.
            Defaults to 1e-4.

    Returns:
        A ListedColormap with transparency applied.

    """
    cmap_with_alpha = cmap(np.arange(cmap.N))

    values = np.linspace(0, 1, cmap.N)
    alpha = 1 - np.log(values[1:]) / np.log(min_value)
    cmap_with_alpha[0, -1] = 0
    cmap_with_alpha[1:, -1] = np.clip(alpha, 0, 1)

    return ListedColormap(cmap_with_alpha)


alphamagma = _colormap_add_transparency(plt.get_cmap("magma_r"), min_value=1e-4)
