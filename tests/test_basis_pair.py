# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

if TYPE_CHECKING:
    from pairinteraction import BasisPair

    from .utils import PairinteractionModule


@pytest.fixture
def basis(pi_module: PairinteractionModule) -> BasisPair:
    """Create a test basis with a few states around Rb 60S."""
    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    energy_min = ket.get_energy(unit="GHz") - 100
    energy_max = ket.get_energy(unit="GHz") + 100
    basis_atom = pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 2), energy=(energy_min, energy_max), energy_unit="GHz")
    system_atom = pi_module.SystemAtom(basis_atom).set_electric_field([0.1, 0, 0.2], unit="V/cm")
    system_atom.diagonalize()
    return pi_module.BasisPair([system_atom, system_atom])


def test_basis_creation(basis: BasisPair) -> None:
    """Test basic properties of created basis."""
    assert basis.number_of_kets == 80 * 80
    assert basis.number_of_states == basis.number_of_kets
    assert len(basis.kets) == basis.number_of_kets
    assert all(x in str(basis) for x in ["BasisPair", "Rb:58,S_1/2,", "...", "Rb:61,D_5/2,"])


def test_coefficients(basis: BasisPair) -> None:
    """Test coefficient matrix properties."""
    coeffs = basis.get_coefficients()
    assert coeffs.shape == (basis.number_of_kets, basis.number_of_states)
    assert pytest.approx(coeffs.diagonal()) == 1.0  # NOSONAR
    assert pytest.approx(coeffs.sum()) == basis.number_of_kets  # NOSONAR


def test_get_amplitudes_and_overlaps(basis: BasisPair) -> None:
    """Test amplitude and overlap calculations."""
    # Test with ket
    test_ket = basis.kets[0]
    amplitudes = basis.get_amplitudes(test_ket)
    assert len(amplitudes) == basis.number_of_states
    assert pytest.approx(amplitudes[0]) == 1.0  # NOSONAR
    overlaps = basis.get_overlaps(test_ket)
    assert len(overlaps) == basis.number_of_states
    assert pytest.approx(overlaps[0]) == 1.0  # NOSONAR

    # Test with basis
    matrix_amplitudes = basis.get_amplitudes(basis)
    assert matrix_amplitudes.shape == (basis.number_of_kets, basis.number_of_states)
    assert pytest.approx(matrix_amplitudes.diagonal()) == 1.0  # NOSONAR
    matrix_overlaps = basis.get_overlaps(basis)
    assert matrix_overlaps.shape == (basis.number_of_states, basis.number_of_states)
    assert pytest.approx(matrix_overlaps.diagonal()) == 1.0  # NOSONAR


def test_get_matrix_elements(basis: BasisPair) -> None:
    """Test matrix element calculations."""
    # Test with ket
    test_ket = basis.kets[0]
    elements_dipole = basis.get_matrix_elements(
        test_ket, ("electric_dipole", "electric_dipole"), qs=(0, 0), unit="e^2 * a0^2"
    )
    assert elements_dipole.shape == (basis.number_of_states,)
    assert np.count_nonzero(elements_dipole) > 0
    assert np.count_nonzero(elements_dipole) <= basis.number_of_states

    # Test with basis
    matrix_elements = basis.get_matrix_elements(
        basis, ("electric_dipole", "electric_dipole"), qs=(0, 0), unit="e^2 * a0^2"
    )
    assert matrix_elements.shape == (basis.number_of_states, basis.number_of_states)
    assert np.count_nonzero(matrix_elements.toarray()) > 0
    assert np.count_nonzero(matrix_elements.toarray()) <= basis.number_of_states**2


def test_error_handling(basis: BasisPair) -> None:
    """Test error cases."""
    with pytest.raises(TypeError):
        basis.get_amplitudes("not a ket")  # type: ignore [arg-type]

    with pytest.raises(TypeError):
        basis.get_overlaps("not a ket")  # type: ignore [arg-type]

    with pytest.raises(TypeError):
        basis.get_matrix_elements("not a ket", ("energy", "energy"), (0, 0))  # type: ignore [arg-type]
