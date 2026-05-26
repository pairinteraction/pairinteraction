# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from pairinteraction import BasisAtom
from pairinteraction.ket.ket_atom import KetAtom
from pairinteraction.state.state_atom import StateAtom
from scipy.sparse import csr_matrix

if TYPE_CHECKING:
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


def test_restriction_mode(pi_module: PairinteractionModule) -> None:
    """Test exact, fuzzy, and numeric restrictions for expectation-value quantum numbers."""
    l_range = (1, 1)
    fuzzy_basis = pi_module.BasisAtom("Yb171_mqdt", nu=(58, 62), l=l_range, m=(0.5, 0.5))
    assert fuzzy_basis.number_of_kets > 0
    assert any(ket.l > l_range[1] for ket in fuzzy_basis.kets)
    assert any(ket.l < l_range[0] for ket in fuzzy_basis.kets)

    exact_basis = pi_module.BasisAtom("Yb171_mqdt", nu=(58, 62), l=l_range, m=(0.5, 0.5), mode="exact")
    assert exact_basis.number_of_kets > 0
    assert all(l_range[0] <= ket.l <= l_range[1] for ket in exact_basis.kets)
    assert exact_basis.number_of_kets < fuzzy_basis.number_of_kets

    factor_basis = pi_module.BasisAtom("Yb171_mqdt", nu=(58, 62), l=l_range, m=(0.5, 0.5), mode=100)
    assert factor_basis.number_of_kets > fuzzy_basis.number_of_kets

    with pytest.raises(ValueError, match="mode"):
        pi_module.BasisAtom("Yb171_mqdt", nu=(58, 62), l=l_range, m=(0.5, 0.5), mode="invalid")  # type: ignore[arg-type]

    with pytest.raises(ValueError, match="non-negative"):
        pi_module.BasisAtom("Yb171_mqdt", nu=(58, 62), l=l_range, m=(0.5, 0.5), mode=-1)


def test_coefficients(basis: BasisAtom) -> None:
    """Test coefficient matrix properties."""
    coeffs = basis.get_coefficients()
    assert coeffs.shape == (basis.number_of_kets, basis.number_of_states)
    assert pytest.approx(coeffs.diagonal()) == 1.0  # NOSONAR
    assert pytest.approx(coeffs.sum()) == basis.number_of_kets  # NOSONAR


def _get_expected_shape(other: KetAtom | StateAtom | BasisAtom, target_basis: BasisAtom) -> tuple[int, ...]:
    if isinstance(other, (KetAtom, StateAtom)):
        return (target_basis.number_of_states,)
    if isinstance(other, BasisAtom):
        return (other.number_of_states, target_basis.number_of_states)
    raise ValueError("Invalid basis_like type")


@pytest.mark.parametrize(
    "other_key",
    [
        "ket_from_basis",
        "ket_from_basis2",
        "other_ket_from_basis2",
        "state_from_basis",
        "state_from_basis2",
        "basis",
        "basis2",
    ],
)
def test_get_methods(basis: BasisAtom, basis2: BasisAtom, other_key: str) -> None:
    """Test amplitude, overlap and matrix element calculations with another ket, state and basis."""
    ind = 5
    other_dict: dict[str, KetAtom | StateAtom | BasisAtom] = {
        "ket_from_basis": basis.get_ket(ind),
        "ket_from_basis2": next(ket for ket in basis2.kets if ket in basis.kets),
        "other_ket_from_basis2": next(ket for ket in basis2.kets if ket not in basis.kets and ket.j_ryd < 3),
        "state_from_basis": basis.get_state(ind),
        "state_from_basis2": basis2.get_state(0),
        "basis": basis,
        "basis2": basis2,
    }
    other = other_dict[other_key]

    amplitudes = basis.get_amplitudes(other)
    assert amplitudes.shape == _get_expected_shape(other, basis)

    overlaps = basis.get_overlaps(other)
    assert overlaps.shape == _get_expected_shape(other, basis)

    matrix_elements = basis.get_matrix_elements(other, "electric_dipole", q=0, unit="e * a0")
    assert matrix_elements.shape == _get_expected_shape(other, basis)

    if other_key in ["ket_from_basis", "state_from_basis"]:
        assert pytest.approx(amplitudes[ind]) == 1.0  # NOSONAR
        assert pytest.approx(overlaps[ind]) == 1.0  # NOSONAR

    if other_key == "ket_from_basis2":
        assert isinstance(amplitudes, np.ndarray)
        assert isinstance(overlaps, np.ndarray)
        assert pytest.approx(np.max(amplitudes)) == 1.0  # NOSONAR
        assert pytest.approx(np.max(overlaps)) == 1.0  # NOSONAR

    if other_key == "other_ket_from_basis2":
        assert isinstance(amplitudes, np.ndarray)
        assert isinstance(overlaps, np.ndarray)
        assert np.count_nonzero(amplitudes) == 0
        assert np.count_nonzero(overlaps) == 0

    if other_key == "basis":
        assert isinstance(amplitudes, csr_matrix)
        assert isinstance(overlaps, csr_matrix)
        assert pytest.approx(amplitudes.diagonal()) == 1.0  # NOSONAR
        assert pytest.approx(overlaps.diagonal()) == 1.0  # NOSONAR

    if other_key == "basis2":
        assert isinstance(amplitudes, csr_matrix)
        assert isinstance(overlaps, csr_matrix)
        n_matching_kets = len([ket for ket in basis2.kets if ket in basis.kets])
        assert np.count_nonzero(amplitudes.toarray()) == n_matching_kets
        assert np.count_nonzero(overlaps.toarray()) == n_matching_kets

    if other_key.startswith("basis"):
        assert isinstance(matrix_elements, csr_matrix)
        assert 0 < np.count_nonzero(matrix_elements.toarray()) < basis.number_of_states**2
    else:
        assert isinstance(matrix_elements, np.ndarray)
        assert 0 < np.count_nonzero(matrix_elements) < basis.number_of_states


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
