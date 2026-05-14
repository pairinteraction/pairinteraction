# SPDX-FileCopyrightText: 2026 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from pairinteraction.units import QuantityScalar, ureg

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import TypeAlias

    from pairinteraction.units import NDArray, PintFloat

    PermittivityLike: TypeAlias = "complex | Callable[[PintFloat], complex]"


def evaluate_relative_permittivity(
    epsilon: PermittivityLike, transition_energy: float, transition_energy_unit: str | None = None
) -> complex:
    """Get the electric permittivity for the given frequency.

    Args:
        epsilon: The electric permittivity (dimensionless) or a callable function that returns it.
        transition_energy: The angular frequency at which to evaluate the permittivity.
            Only needed if the permittivity is frequency dependent.
        transition_energy_unit: The unit of the angular frequency.
            Default None, which means that the angular frequency must be given as pint object.

    Returns:
        The electric permittivity at the given angular frequency.

    """
    omega_au = QuantityScalar.convert_user_to_au(transition_energy, transition_energy_unit, "energy")
    if np.isscalar(epsilon):
        return epsilon  # type: ignore [return-value]
    if callable(epsilon):
        return epsilon(ureg.Quantity(omega_au, "hartree"))
    raise TypeError("epsilon must be either a complex number or a callable function.")


def normalize(vector: NDArray) -> NDArray:
    """Return the normalized version of the input vector."""
    norm = np.linalg.norm(vector)
    if np.isclose(norm, 0):
        raise ValueError("Cannot normalize a zero vector.")
    return vector / norm  # type: ignore [no-any-return]


def get_lab_to_local_rotation_matrix(normal: NDArray) -> NDArray:
    """Return a rotation matrix that maps lab coordinates to a frame with local z parallel to normal."""
    if np.isclose(np.linalg.norm(normal), 0):
        raise ValueError("Normal vector cannot be zero.")

    z_axis = normalize(normal)
    reference_axis = np.eye(3)[np.argmin(np.abs(z_axis))]
    x_axis = normalize(np.cross(reference_axis, z_axis))
    y_axis = normalize(np.cross(z_axis, x_axis))

    return np.vstack((x_axis, y_axis, z_axis))


def rotate_vector_to_local(vector_lab: NDArray, lab_to_local_rotation: NDArray) -> NDArray:
    return lab_to_local_rotation @ vector_lab


def rotate_tensor_to_lab(tensor_local: NDArray, lab_to_local_rotation: NDArray) -> NDArray:
    return lab_to_local_rotation.T @ tensor_local @ lab_to_local_rotation
