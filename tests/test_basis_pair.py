# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

if TYPE_CHECKING:
    from pairinteraction import BasisAtom, BasisPair, KetAtom, SystemAtom

    from .utils import PairinteractionModule

from .utils import no_log_propagation


@pytest.fixture
def basis_atom(pi_module: PairinteractionModule) -> BasisAtom:
    """Create a BasisAtom around Rb 60S."""
    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    energy_min = ket.get_energy(unit="GHz") - 100
    energy_max = ket.get_energy(unit="GHz") + 100
    return pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 2), energy=(energy_min, energy_max), energy_unit="GHz")


@pytest.fixture
def system_atom(pi_module: PairinteractionModule, basis_atom: BasisAtom) -> SystemAtom:
    """Create a diagonalized SystemAtom around Rb 60S."""
    system_atom = pi_module.SystemAtom(basis_atom)
    system_atom.diagonalize()
    return system_atom


@pytest.fixture
def basis(pi_module: PairinteractionModule, system_atom: SystemAtom) -> BasisPair:
    """Create a test basis with a few states around Rb 60S."""
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


def test_get_amplitudes_and_overlaps(pi_module: PairinteractionModule, basis: BasisPair) -> None:
    """Test amplitude and overlap calculations."""
    # Test with ket
    test_ket = basis.get_ket(0)
    amplitudes = basis.get_amplitudes(test_ket)
    assert len(amplitudes) == basis.number_of_states
    assert pytest.approx(amplitudes[0]) == 1.0  # NOSONAR
    overlaps = basis.get_overlaps(test_ket)
    assert len(overlaps) == basis.number_of_states
    assert pytest.approx(overlaps[0]) == 1.0  # NOSONAR

    # Test with ket atom tuple
    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    overlaps = basis.get_overlaps((ket, ket))
    assert overlaps.shape == (basis.number_of_states,)

    # Test with basis
    matrix_amplitudes = basis.get_amplitudes(basis)
    assert matrix_amplitudes.shape == (basis.number_of_kets, basis.number_of_states)
    assert pytest.approx(matrix_amplitudes.diagonal()) == 1.0  # NOSONAR
    matrix_overlaps = basis.get_overlaps(basis)
    assert matrix_overlaps.shape == (basis.number_of_states, basis.number_of_states)
    assert pytest.approx(matrix_overlaps.diagonal()) == 1.0  # NOSONAR

    # Test with a different basis (non-identical kets) - cross-basis overlaps
    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    energy_min = ket.get_energy(unit="GHz") - 22
    energy_max = ket.get_energy(unit="GHz") + 22
    basis_atom_new = pi_module.BasisAtom("Rb", n=(60, 62), l=(0, 1), energy=(energy_min, energy_max), energy_unit="GHz")
    system_atom_new = pi_module.SystemAtom(basis_atom_new)
    basis2 = pi_module.BasisPair([system_atom_new, system_atom_new])

    matrix_overlaps = basis.get_overlaps(basis2)
    assert matrix_overlaps.shape == (basis2.number_of_states, basis.number_of_states)
    assert (matrix_overlaps.data >= 0).all()
    assert (matrix_overlaps.data <= 1 + 1e-10).all()

    col_sums = np.array(matrix_overlaps.sum(axis=0)).flatten()
    assert np.all(col_sums <= 1.0 + 1e-10)
    row_sums = np.array(matrix_overlaps.sum(axis=1)).flatten()
    assert np.all(row_sums <= 1.0 + 1e-10)

    idx0 = basis.get_corresponding_state_index((ket, ket))
    idx2 = basis2.get_corresponding_state_index((ket, ket))
    assert pytest.approx(matrix_overlaps[idx2, idx0]) == 1.0  # NOSONAR


def test_get_matrix_elements(basis: BasisPair) -> None:
    """Test matrix element calculations."""
    # Test with ket
    test_ket = basis.get_ket(0)
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


def test_from_ket_atoms(pi_module: PairinteractionModule, system_atom: SystemAtom) -> None:
    """Test BasisPair.from_ket_atoms."""
    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket1 = pi_module.KetAtom("Rb", n=59, l=1, j=0.5, m=0.5)
    ket2 = pi_module.KetAtom("Rb", n=60, l=1, j=0.5, m=0.5)
    ket_atoms_dict: dict[str, tuple[KetAtom, KetAtom] | list[tuple[KetAtom, KetAtom]]] = {
        "single_ket": (ket, ket),
        "multiple_kets": [(ket, ket), (ket1, ket2), (ket2, ket1)],
    }

    for ket_atoms in ket_atoms_dict.values():
        # delta_energy restriction
        pair_basis = pi_module.BasisPair.from_ket_atoms(
            ket_atoms, [system_atom, system_atom], delta_energy=3, delta_energy_unit="GHz", delta_m=1
        )
        assert pair_basis.number_of_kets > 0

        # number_of_kets
        ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
        for target in [50, 100, 1000]:
            with no_log_propagation("pairinteraction.basis.basis_pair"):  # surpress number_of_kets warning
                pair_basis = pi_module.BasisPair.from_ket_atoms(
                    ket_atoms, [system_atom, system_atom], number_of_kets=target, delta_m=1
                )
            assert pair_basis.number_of_kets >= target
            assert pair_basis.number_of_kets < target + 20  # allow some extra due to degeneracies

    # test error cases
    with pytest.raises(ValueError, match="empty"):
        pi_module.BasisPair.from_ket_atoms([], [system_atom, system_atom])

    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    with pytest.raises(ValueError, match="number_of_kets"):
        pi_module.BasisPair.from_ket_atoms(
            (ket, ket),
            [system_atom, system_atom],
            delta_energy=3,
            delta_energy_unit="GHz",
            number_of_kets=10,
        )
