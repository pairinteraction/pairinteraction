# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from pairinteraction import BasisAtom, BasisPair
from pairinteraction.ket.ket_pair import is_ket_pair_like
from pairinteraction.state.state_pair import is_state_pair_like

if TYPE_CHECKING:
    from pairinteraction import KetAtom, SystemAtom
    from pairinteraction.basis.basis_pair import BasisPairLike
    from pairinteraction.ket.ket_pair import KetPairLike
    from pairinteraction.state.state_pair import StatePairLike

    from .utils import PairinteractionModule

from .utils import no_log_propagation


@pytest.fixture
def basis(pi_module: PairinteractionModule, system_atom: SystemAtom, system_atom2: SystemAtom) -> BasisPair:
    system_atoms = [system_atom, system_atom2]
    ket_atoms = tuple(s.basis.get_ket(30) for s in system_atoms)
    with no_log_propagation("pairinteraction.basis.basis_pair"):  # suppress number_of_kets warning
        return pi_module.BasisPair.from_ket_atoms(ket_atoms, system_atoms, number_of_kets=80)


@pytest.fixture
def basis2(pi_module: PairinteractionModule, system_atom: SystemAtom) -> BasisPair:
    system_atoms = [system_atom, system_atom]
    ket_atoms = tuple(s.basis.get_ket(30) for s in system_atoms)
    with no_log_propagation("pairinteraction.basis.basis_pair"):  # suppress number_of_kets warning
        return pi_module.BasisPair.from_ket_atoms(ket_atoms, system_atoms, number_of_kets=80)


@pytest.fixture
def basis_atom(pi_module: PairinteractionModule) -> BasisAtom:
    return pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 2))


@pytest.fixture
def basis_atom2(pi_module: PairinteractionModule) -> BasisAtom:
    return pi_module.BasisAtom("Rb", n=(58, 62), l=(2, 3))


@pytest.fixture
def system_atom(pi_module: PairinteractionModule, basis_atom: BasisAtom) -> SystemAtom:
    return pi_module.SystemAtom(basis_atom).diagonalize()


@pytest.fixture
def system_atom2(pi_module: PairinteractionModule, basis_atom2: BasisAtom) -> SystemAtom:
    return pi_module.SystemAtom(basis_atom2).diagonalize()


def test_basis_creation(pi_module: PairinteractionModule, system_atom: SystemAtom, system_atom2: SystemAtom) -> None:
    """Test basic properties of created basis."""
    basis = pi_module.BasisPair([system_atom, system_atom2])
    assert basis.number_of_kets == system_atom.basis.number_of_kets * system_atom2.basis.number_of_kets
    assert basis.number_of_states == basis.number_of_kets
    assert len(basis.kets) == basis.number_of_kets
    assert all(x in str(basis) for x in ["BasisPair", "Rb:58,S_1/2,", "...", "Rb:62,D_5/2,"])


def test_coefficients(basis: BasisPair) -> None:
    """Test coefficient matrix properties."""
    coeffs = basis.get_coefficients()
    assert coeffs.shape == (basis.number_of_kets, basis.number_of_states)
    assert pytest.approx(coeffs.diagonal()) == 1.0  # NOSONAR
    assert pytest.approx(coeffs.sum()) == basis.number_of_kets  # NOSONAR


def _get_expected_shape(other: KetPairLike | StatePairLike | BasisPairLike, target_basis: BasisPair) -> tuple[int, ...]:
    if is_ket_pair_like(other) or is_state_pair_like(other):
        return (target_basis.number_of_states,)
    if isinstance(other, BasisPair):
        return (other.number_of_states, target_basis.number_of_states)
    if isinstance(other, list) and all(isinstance(b, BasisAtom) for b in other):
        return (other[0].number_of_states * other[1].number_of_states, target_basis.number_of_states)
    raise ValueError("Invalid basis_like type")


@pytest.mark.parametrize(
    "other_key",
    [
        "ket_from_basis",
        "ket_from_basis2",
        "ket_atom_tuple_from_basis",
        "ket_atom_tuple_from_basis2",
        "state_from_basis",
        "state_from_basis2",
        "state_atom_tuple_from_basis",
        "state_atom_tuple_from_basis2",
        "basis",
        "basis2",
        "basis_atom_tuple_from_basis",
        "basis_atom_tuple_from_basis2",
    ],
)
def test_get_methods(basis: BasisPair, basis2: BasisPair, other_key: str) -> None:
    """Test amplitude, overlap and matrix element calculations with another ket, state and basis."""
    other_dict: dict[str, KetPairLike | StatePairLike | BasisPairLike] = {
        "ket_from_basis": basis.get_ket(0),
        "ket_from_basis2": basis2.get_ket(0),
        "ket_atom_tuple_from_basis": [s.basis.get_ket(0) for s in basis.system_atoms],
        "ket_atom_tuple_from_basis2": [s.basis.get_ket(0) for s in basis2.system_atoms],
        "state_from_basis": basis.get_state(0),
        "state_from_basis2": basis2.get_state(0),
        "state_atom_tuple_from_basis": [s.basis.get_state(0) for s in basis.system_atoms],
        "state_atom_tuple_from_basis2": [s.basis.get_state(0) for s in basis2.system_atoms],
        "basis": basis,
        "basis2": basis2,
        "basis_atom_tuple_from_basis": [s.basis for s in basis.system_atoms],
        "basis_atom_tuple_from_basis2": [s.basis for s in basis2.system_atoms],
    }
    other = other_dict[other_key]

    amplitudes = basis.get_amplitudes(other)
    assert amplitudes.shape == _get_expected_shape(other, basis)

    overlaps = basis.get_overlaps(other)
    assert overlaps.shape == _get_expected_shape(other, basis)

    elements_dipole = basis.get_matrix_elements(
        other, ("electric_dipole", "electric_dipole"), qs=(0, 0), unit="e^2 * a0^2"
    )
    assert elements_dipole.shape == _get_expected_shape(other, basis)


def test_get_overlaps_pair_explicitly(pi_module: PairinteractionModule, system_atom: SystemAtom) -> None:
    # Test with a different basis (non-identical kets) - cross-basis overlaps
    ket = system_atom.basis.get_ket(20)
    energy_min = 2 * ket.get_energy("GHz") - 1
    energy_max = 2 * ket.get_energy("GHz") + 1
    basis_pair1 = pi_module.BasisPair([system_atom, system_atom])
    basis_pair2 = pi_module.BasisPair([system_atom, system_atom], energy=(energy_min, energy_max), energy_unit="GHz")
    assert basis_pair1.number_of_states > basis_pair2.number_of_states
    assert basis_pair1.number_of_states > 1

    matrix_overlaps = basis_pair1.get_overlaps(basis_pair2)
    assert matrix_overlaps.shape == (basis_pair2.number_of_states, basis_pair1.number_of_states)
    assert (matrix_overlaps.data >= 0).all()
    assert (matrix_overlaps.data <= 1 + 1e-10).all()

    col_sums = np.array(matrix_overlaps.sum(axis=0)).flatten()
    assert np.all(col_sums <= 1.0 + 1e-10)
    row_sums = np.array(matrix_overlaps.sum(axis=1)).flatten()
    assert np.all(row_sums <= 1.0 + 1e-10)

    idx0 = basis_pair1.get_corresponding_state_index((ket, ket))
    idx2 = basis_pair2.get_corresponding_state_index((ket, ket))
    assert pytest.approx(matrix_overlaps[idx2, idx0]) == 1.0  # NOSONAR


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
            with no_log_propagation("pairinteraction.basis.basis_pair"):  # suppress number_of_kets warning
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
