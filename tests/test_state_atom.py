# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

if TYPE_CHECKING:
    from pairinteraction import BasisAtom, StateAtom

    from .utils import PairinteractionModule


@pytest.fixture
def basis(pi_module: PairinteractionModule) -> BasisAtom:
    """Create a test basis with a few states around Rb 60S."""
    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    energy_min = ket.get_energy(unit="GHz") - 100
    energy_max = ket.get_energy(unit="GHz") + 100
    return pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 2), energy=(energy_min, energy_max), energy_unit="GHz")


@pytest.fixture
def state(pi_module: PairinteractionModule, basis: BasisAtom) -> StateAtom:
    """Create a test state."""
    ket = pi_module.KetAtom("Rb", n=60, l=1, j=1.5, m=-0.5)
    return basis.get_corresponding_state(ket)


def test_state_creation(state: StateAtom) -> None:
    """Test basic properties of created state."""
    assert state.species == "Rb"
    assert state.number_of_kets == 80
    assert len(state.kets) == state.number_of_kets
    assert all(x in str(state) for x in ["StateAtom", "60", "S", "3/2", "-1/2"])
    assert state.is_canonical


def test_coefficients(state: StateAtom) -> None:
    """Test coefficient matrix properties."""
    coeffs = state.get_coefficients()
    assert coeffs.shape == (state.number_of_kets,)
    assert np.count_nonzero(coeffs) == 1
    assert pytest.approx(coeffs.sum()) == 1.0  # NOSONAR


def test_get_amplitude_and_overlap(state: StateAtom) -> None:
    """Test amplitude and overlap calculations."""
    # Test with ket
    test_ket = state.get_corresponding_ket()
    amplitude = state.get_amplitude(test_ket)
    assert np.isscalar(amplitude)
    assert pytest.approx(amplitude) == 1.0  # NOSONAR
    overlap = state.get_overlap(test_ket)
    assert np.isscalar(overlap)
    assert pytest.approx(overlap) == 1.0  # NOSONAR

    # Test with state
    amplitude = state.get_amplitude(state)
    assert np.isscalar(amplitude)
    assert pytest.approx(amplitude) == 1.0  # NOSONAR
    overlap = state.get_overlap(state)
    assert np.isscalar(overlap)
    assert pytest.approx(overlap) == 1.0  # NOSONAR


def test_get_matrix_element(pi_module: PairinteractionModule, basis: BasisAtom) -> None:
    """Test matrix element calculations."""
    ket1 = pi_module.KetAtom("Rb", n=60, l=1, j=1.5, m=-0.5)
    ket2 = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    state1 = basis.get_corresponding_state(ket1)
    state2 = basis.get_corresponding_state(ket2)

    # Test with ket
    element_dipole_ket = state1.get_matrix_element(ket2, "electric_dipole", q=1, unit="e * a0")
    assert np.isscalar(element_dipole_ket)
    assert element_dipole_ket != 0

    # Test with state
    element_dipole_state = state1.get_matrix_element(state2, "electric_dipole", q=1, unit="e * a0")
    assert np.isscalar(element_dipole_state)
    assert pytest.approx(element_dipole_ket) == element_dipole_state  # NOSONAR
    assert state1.get_matrix_element(state1, "electric_dipole", q=1, unit="e * a0") == 0
    assert state1.get_matrix_element(state2, "electric_dipole", q=0, unit="e * a0") == 0


def test_state_without_basis(pi_module: PairinteractionModule) -> None:
    """A StateAtom created without an explicit basis uses a minimal single-ket basis."""
    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    state = pi_module.StateAtom(ket)
    assert state.number_of_kets == 1
    assert state.is_canonical
    assert state.get_corresponding_ket() == ket
    assert pytest.approx(state.get_overlap(ket)) == 1.0  # NOSONAR


def test_add_states_with_same_basis(basis: BasisAtom, pi_module: PairinteractionModule) -> None:
    """Adding states expressed in the same basis keeps that basis."""
    ket1 = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket2 = pi_module.KetAtom("Rb", n=60, l=1, j=1.5, m=0.5)
    state1 = basis.get_corresponding_state(ket1)
    state2 = basis.get_corresponding_state(ket2)

    combined = (state1 + state2).normalize()
    assert combined.number_of_kets == basis.number_of_kets
    assert pytest.approx(combined.get_overlap(ket1)) == 0.5  # NOSONAR
    assert pytest.approx(combined.get_overlap(ket2)) == 0.5  # NOSONAR


def test_add_states_with_different_bases_merges(pi_module: PairinteractionModule) -> None:
    """Adding states with different bases merges the bases and re-expresses the coefficients."""
    ket1 = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket2 = pi_module.KetAtom("Rb", n=60, l=1, j=1.5, m=0.5)

    # Both states start from a minimal single-ket basis, so their bases differ.
    state1 = pi_module.StateAtom(ket1)
    state2 = pi_module.StateAtom(ket2)
    assert state1.number_of_kets == 1
    assert state2.number_of_kets == 1

    combined = (state1 + 2 * state2).normalize()
    # The merged basis contains both kets.
    assert combined.number_of_kets == 2
    assert pytest.approx(combined.get_overlap(ket1)) == 1 / 5  # NOSONAR
    assert pytest.approx(combined.get_overlap(ket2)) == 4 / 5  # NOSONAR

    # Adding a single-ket state onto a full-basis state merges into the full basis.
    energy_min = ket1.get_energy(unit="GHz") - 100
    energy_max = ket1.get_energy(unit="GHz") + 100
    full_basis = pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 2), energy=(energy_min, energy_max), energy_unit="GHz")
    state_full = full_basis.get_corresponding_state(ket1)
    merged = (state_full + pi_module.StateAtom(ket2)).normalize()
    assert merged.number_of_kets == full_basis.number_of_kets
    assert pytest.approx(merged.get_overlap(ket1)) == 0.5  # NOSONAR
    assert pytest.approx(merged.get_overlap(ket2)) == 0.5  # NOSONAR


def test_error_handling(state: StateAtom) -> None:
    """Test error cases."""
    with pytest.raises(TypeError):
        state.get_amplitude("not a ket")  # type: ignore [arg-type]

    with pytest.raises(TypeError):
        state.get_overlap("not a ket")  # type: ignore [arg-type]

    with pytest.raises(TypeError):
        state.get_matrix_element("not a ket", "energy", 0)  # type: ignore [call-overload]
