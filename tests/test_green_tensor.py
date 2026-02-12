# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from pairinteraction import ureg
from pairinteraction.green_tensor import GreenTensorFreeSpace
from pairinteraction.green_tensor.green_tensor_cavity import GreenTensorCavity
from pairinteraction.green_tensor.green_tensor_surface import GreenTensorSurface

if TYPE_CHECKING:
    from pairinteraction.green_tensor.green_tensor_base import GreenTensorBase
    from pairinteraction.units import NDArray

    from .utils import PairinteractionModule

DISTANCE_VECTOR_MUM_LIST = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 3.5]]


@pytest.mark.parametrize("distance_vector_mum", DISTANCE_VECTOR_MUM_LIST)
def test_static_green_tensor(distance_vector_mum: list[float]) -> None:
    gt_reference = reference_green_tensor_vacuum(distance_vector_mum)

    gt: GreenTensorBase
    gt = GreenTensorFreeSpace([0, 0, 0], distance_vector_mum, unit="micrometer", static_limit=True)
    gt_dipole_dipole = gt.get(1, 1, transition_energy=0, scaled=True)
    np.testing.assert_allclose(gt_dipole_dipole.m, gt_reference)

    gt = GreenTensorSurface([0, 0, 0], distance_vector_mum, z=1000, unit="micrometer", static_limit=True)
    gt_dipole_dipole = gt.get(1, 1, transition_energy=0, scaled=True)
    np.testing.assert_allclose(gt_dipole_dipole.m, gt_reference, rtol=1e-6, atol=1e-16)

    gt = GreenTensorCavity([0, 0, 0], distance_vector_mum, z1=-1000, z2=1000, unit="micrometer", static_limit=True)
    gt_dipole_dipole = gt.get(1, 1, transition_energy=0, scaled=True)
    np.testing.assert_allclose(gt_dipole_dipole.m, gt_reference, rtol=1e-6, atol=1e-16)


@pytest.mark.parametrize("distance_vector_mum", DISTANCE_VECTOR_MUM_LIST)
def test_static_green_tensor_pair_potential(pi_module: PairinteractionModule, distance_vector_mum: list[float]) -> None:
    ket = pi_module.KetAtom("Rb", n=60, l=0, m=0.5)
    basis = pi_module.BasisAtom("Rb", n=(ket.n - 2, ket.n + 2), l=(0, 2))
    system = pi_module.SystemAtom(basis)
    system.set_magnetic_field([0, 0, 1], "gauss")
    system.set_electric_field([0, 0, 1], "V/cm")
    system.diagonalize()

    pair_energy_ghz = 2 * system.get_corresponding_energy(ket, unit="GHz")
    energy_range_ghz = (pair_energy_ghz - 3, pair_energy_ghz + 3)
    basis_pair = pi_module.BasisPair([system, system], energy=energy_range_ghz, energy_unit="GHz")

    system_pair_vacuum = pi_module.SystemPair(basis_pair)
    system_pair_vacuum.set_distance_vector(distance_vector_mum, unit="micrometer")
    system_pair_vacuum.set_interaction_order(3)

    gt = GreenTensorFreeSpace([0, 0, 0], distance_vector_mum, unit="micrometer", interaction_order=3, static_limit=True)
    gt.set_relative_permittivity(1.0)
    system_pair_free_space = pi_module.SystemPair(basis_pair).set_green_tensor(gt)

    pi_module.diagonalize([system_pair_free_space, system_pair_vacuum], sort_by_energy=True)
    np.testing.assert_allclose(
        system_pair_vacuum.get_eigenenergies("GHz") - pair_energy_ghz,
        system_pair_free_space.get_eigenenergies("GHz") - pair_energy_ghz,
    )


@pytest.mark.parametrize("distance_vector_mum", DISTANCE_VECTOR_MUM_LIST)
def test_vacuum_green_tensor(pi_module: PairinteractionModule, distance_vector_mum: list[float]) -> None:
    ket1 = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket2 = pi_module.KetAtom("Rb", n=60, l=1, j=0.5, m=0.5)

    basis = pi_module.BasisAtom("Rb", n=(0, 0), additional_kets=[ket1, ket2])
    system = pi_module.SystemAtom(basis)
    pair_energy = ket1.get_energy("GHz") + ket2.get_energy("GHz")
    basis_pair = pi_module.BasisPair((system, system), energy=(pair_energy - 0.1, pair_energy + 0.1), energy_unit="GHz")

    dd = ket1.get_matrix_element(ket2, "electric_dipole", q=0)
    gt_reference = reference_green_tensor_vacuum(distance_vector_mum)
    reference = dd * gt_reference[2, 2] * dd

    # test internal vacuum green tensor
    system_pair = pi_module.SystemPair(basis_pair)
    system_pair.set_distance_vector(distance_vector_mum, "micrometer")
    hamiltonian = system_pair.get_hamiltonian("hartree").toarray()
    np.testing.assert_allclose(hamiltonian[1, 0], reference)

    # test custom free space vacuum green tensor
    system_pair = pi_module.SystemPair(basis_pair)
    gt = GreenTensorFreeSpace([0, 0, 0], distance_vector_mum, unit="micrometer", static_limit=True)
    system_pair.set_green_tensor(gt)
    hamiltonian = system_pair.get_hamiltonian("hartree").toarray()
    np.testing.assert_allclose(hamiltonian[1, 0], reference)


def reference_green_tensor_vacuum(distance_vector_mum: list[float]) -> NDArray:
    distance_mum = np.linalg.norm(distance_vector_mum)
    distance_au = ureg.Quantity(distance_mum, "micrometer").to_base_units().m
    return (np.eye(3) - 3 * np.outer(distance_vector_mum, distance_vector_mum) / distance_mum**2) / distance_au**3
