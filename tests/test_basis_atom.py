# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

if TYPE_CHECKING:
    from pairinteraction import BasisAtom

    from .utils import PairinteractionModule


@pytest.fixture
def basis(pi_module: PairinteractionModule) -> BasisAtom:
    return pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 2))


@pytest.fixture
def basis2(pi_module: PairinteractionModule) -> BasisAtom:
    return pi_module.BasisAtom("Rb", n=(58, 62), l=(2, 3))


def test_basis_creation(pi_module: PairinteractionModule) -> None:
    """Test basic properties of created basis."""
    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    energy_min = ket.get_energy(unit="GHz") - 100
    energy_max = ket.get_energy(unit="GHz") + 100
    basis = pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 2), energy=(energy_min, energy_max), energy_unit="GHz")
    assert basis.species == "Rb"
    assert basis.number_of_kets == 80
    assert basis.number_of_states == basis.number_of_kets
    assert len(basis.kets) == basis.number_of_kets
    assert basis.number_of_kets < pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 2)).number_of_kets
    assert all(x in str(basis) for x in ["BasisAtom", "n=(58, 62)", "l=(0, 2)"])


def test_coefficients(basis: BasisAtom) -> None:
    """Test coefficient matrix properties."""
    coeffs = basis.get_coefficients()
    assert coeffs.shape == (basis.number_of_kets, basis.number_of_states)
    assert pytest.approx(coeffs.diagonal()) == 1.0  # NOSONAR
    assert pytest.approx(coeffs.sum()) == basis.number_of_kets  # NOSONAR


def test_get_methods_with_ket(basis: BasisAtom, basis2: BasisAtom) -> None:
    """Test amplitude, overlap and matrix element calculations with a ket."""
    # ket from basis
    ket = basis.get_ket(0)
    amplitudes = basis.get_amplitudes(ket)
    assert len(amplitudes) == basis.number_of_states
    assert pytest.approx(amplitudes[0]) == 1.0  # NOSONAR

    overlaps = basis.get_overlaps(ket)
    assert len(overlaps) == basis.number_of_states
    assert pytest.approx(overlaps[0]) == 1.0  # NOSONAR

    elements_dipole = basis.get_matrix_elements(ket, "electric_dipole", q=0, unit="e * a0")
    assert elements_dipole.shape == (basis.number_of_states,)
    assert np.count_nonzero(elements_dipole) > 0
    assert np.count_nonzero(elements_dipole) < basis.number_of_states

    # ket not in basis
    ket = next(basis2.get_ket(i) for i in range(basis2.number_of_kets) if basis2.get_ket(i) not in basis.kets)
    amplitudes = basis.get_amplitudes(ket)
    assert len(amplitudes) == basis.number_of_states
    assert pytest.approx(amplitudes[0]) == 0.0  # NOSONAR

    overlaps = basis.get_overlaps(ket)
    assert len(overlaps) == basis.number_of_states
    assert pytest.approx(overlaps[0]) == 0.0  # NOSONAR

    elements_dipole = basis.get_matrix_elements(ket, "electric_dipole", q=0, unit="e * a0")
    assert elements_dipole.shape == (basis.number_of_states,)
    assert np.count_nonzero(elements_dipole) > 0
    assert np.count_nonzero(elements_dipole) < basis.number_of_states


def test_get_methods_with_state(basis: BasisAtom, basis2: BasisAtom) -> None:
    """Test amplitude, overlap and matrix element calculations with a state."""
    # state from basis
    state = basis.get_state(0)
    amplitudes = basis.get_amplitudes(state)
    assert len(amplitudes) == basis.number_of_states
    assert pytest.approx(amplitudes[0]) == 1.0  # NOSONAR

    overlaps = basis.get_overlaps(state)
    assert len(overlaps) == basis.number_of_states
    assert pytest.approx(overlaps[0]) == 1.0  # NOSONAR

    elements_dipole = basis.get_matrix_elements(state, "electric_dipole", q=0, unit="e * a0")
    assert elements_dipole.shape == (basis.number_of_states,)
    assert np.count_nonzero(elements_dipole) > 0
    assert np.count_nonzero(elements_dipole) < basis.number_of_states

    # state with different basis
    state = basis2.get_state(0)
    amplitudes = basis.get_amplitudes(state)
    assert len(amplitudes) == basis.number_of_states

    overlaps = basis.get_overlaps(state)
    assert len(overlaps) == basis.number_of_states

    elements_dipole = basis.get_matrix_elements(state, "electric_dipole", q=0, unit="e * a0")
    assert elements_dipole.shape == (basis.number_of_states,)


def test_get_methods_with_basis(basis: BasisAtom, basis2: BasisAtom) -> None:
    """Test amplitude, overlap and matrix element calculations with a basis."""
    # Test with same basis
    matrix_amplitudes = basis.get_amplitudes(basis)
    assert matrix_amplitudes.shape == (basis.number_of_kets, basis.number_of_states)
    assert pytest.approx(matrix_amplitudes.diagonal()) == 1.0  # NOSONAR

    matrix_overlaps = basis.get_overlaps(basis)
    assert matrix_overlaps.shape == (basis.number_of_states, basis.number_of_states)
    assert pytest.approx(matrix_overlaps.diagonal()) == 1.0  # NOSONAR

    matrix_elements = basis.get_matrix_elements(basis, "electric_dipole", q=0, unit="e * a0")
    assert matrix_elements.shape == (basis.number_of_states, basis.number_of_states)
    assert np.count_nonzero(matrix_elements.toarray()) > 0
    assert np.count_nonzero(matrix_elements.toarray()) < basis.number_of_states**2

    # Test with other basis
    matrix_amplitudes = basis.get_amplitudes(basis2)
    assert matrix_amplitudes.shape == (basis2.number_of_kets, basis.number_of_states)
    assert pytest.approx(max(matrix_amplitudes.data)) == 1.0  # NOSONAR

    matrix_overlaps = basis.get_overlaps(basis2)
    assert matrix_overlaps.shape == (basis2.number_of_states, basis.number_of_states)
    assert pytest.approx(max(matrix_overlaps.data)) == 1.0  # NOSONAR

    matrix_elements = basis.get_matrix_elements(basis2, "electric_dipole", q=0, unit="e * a0")
    assert matrix_elements.shape == (basis2.number_of_states, basis.number_of_states)
    assert np.count_nonzero(matrix_elements.toarray()) > 0
    assert np.count_nonzero(matrix_elements.toarray()) < basis.number_of_states**2


def test_error_handling(basis: BasisAtom) -> None:
    """Test error cases."""
    with pytest.raises(TypeError):
        basis.get_amplitudes("not a ket")  # type: ignore [call-overload]

    with pytest.raises(TypeError):
        basis.get_overlaps("not a ket")  # type: ignore [call-overload]

    with pytest.raises(TypeError):
        basis.get_matrix_elements("not a ket", "energy", 0)  # type: ignore [call-overload]


def test_from_kets(pi_module: PairinteractionModule) -> None:
    """Test BasisAtom.from_kets."""
    # single ket
    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    basis = pi_module.BasisAtom.from_kets(
        ket,
        delta_n=2,
        delta_nu=3,
        delta_nui=3,
        delta_l=2,
        delta_s=1,
        delta_j=3,
        delta_l_ryd=2,
        delta_j_ryd=3,
        delta_f=3,
        delta_m=2,
        delta_energy=100,
        delta_energy_unit="GHz",
    )
    assert basis.species == "Rb"
    assert all(58 <= k.n <= 62 for k in basis.kets)
    assert any(k.n == 62 for k in basis.kets)
    assert any(k.n == 58 for k in basis.kets)
    assert any(k == ket for k in basis.kets)

    # multiple kets
    ket1 = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket2 = pi_module.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    basis = pi_module.BasisAtom.from_kets([ket1, ket2], delta_n=2)
    assert all(58 <= k.n <= 63 for k in basis.kets)
    assert any(k.n == 63 for k in basis.kets)
    assert any(k.n == 58 for k in basis.kets)
    assert any(k == ket1 for k in basis.kets)
    assert any(k == ket2 for k in basis.kets)

    # test that from_kets is consistent with direct constructor
    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    basis_from = pi_module.BasisAtom.from_kets(ket, delta_n=2)
    basis_direct = pi_module.BasisAtom("Rb", n=(58, 62))
    assert basis_from.number_of_kets == basis_direct.number_of_kets

    # test error cases
    with pytest.raises(ValueError, match="empty"):
        pi_module.BasisAtom.from_kets([])

    ket_rb = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket_sr = pi_module.KetAtom("Sr88_singlet", n=60, l=1, j=1, m=0)
    with pytest.raises(ValueError, match="species"):
        pi_module.BasisAtom.from_kets([ket_rb, ket_sr])
