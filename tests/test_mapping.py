# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Test the mapping between kets and states."""

import numpy as np
import pairinteraction.real as pi
from scipy.optimize import linear_sum_assignment


def test_mapping() -> None:
    """Test generation of a mapping."""
    # Get the eigenbasis of the Hamiltonian describing an atom in an electric field
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2))
    system = pi.SystemAtom(basis).set_electric_field([0, 0, 2.5], unit="V/cm")
    system.diagonalize(diagonalizer="eigen", sort_by_energy=True)
    eigenbasis = system.get_eigenbasis()

    assert eigenbasis.number_of_states == eigenbasis.number_of_kets

    # Obtain the mapping
    state_indices = [eigenbasis.get_corresponding_state_index(ket) for ket in eigenbasis.kets]

    # Calculate the mapping from the coefficient matrix using scipy
    coefficient_matrix = np.square(np.abs(eigenbasis.get_coefficients().todense()))
    rows, cols = linear_sum_assignment(-coefficient_matrix)

    sorter = np.argsort(rows)
    rows = rows[sorter]
    cols = cols[sorter]

    # Because we have chosen the electric field to be weak enough to avoid strong mixing of states,
    # the mapping obtained by pairinteraction's heuristic should be the same as the optimal mapping
    # obtained by scipy's linear_sum_assignment
    np.testing.assert_array_equal(rows, np.arange(eigenbasis.number_of_kets))
    np.testing.assert_array_equal(cols, state_indices)  # TODO
