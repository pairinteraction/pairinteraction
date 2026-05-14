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

    gt = GreenTensorSurface(
        [0, 0, 0],
        distance_vector_mum,
        point_on_plane=[0, 0, 1000],
        surface_normal=[0, 0, 1],
        unit="micrometer",
        static_limit=True,
    )
    gt_dipole_dipole = gt.get(1, 1, transition_energy=0, scaled=True)
    np.testing.assert_allclose(gt_dipole_dipole.m, gt_reference, rtol=1e-6, atol=1e-16)

    gt = GreenTensorCavity(
        [0, 0, 0],
        distance_vector_mum,
        point_on_plane1=[0, 0, -1000],
        point_on_plane2=[0, 0, 1000],
        surface_normal=[0, 0, 1],
        unit="micrometer",
        static_limit=True,
    )
    gt_dipole_dipole = gt.get(1, 1, transition_energy=0, scaled=True)
    np.testing.assert_allclose(gt_dipole_dipole.m, gt_reference, rtol=1e-6, atol=1e-16)


@pytest.mark.parametrize("distance_vector_mum", DISTANCE_VECTOR_MUM_LIST)
def test_static_green_tensor_pair_potential(
    pi_module: PairinteractionModule, use_real: bool, distance_vector_mum: list[float]
) -> None:
    if use_real and distance_vector_mum[1] != 0:
        return  # run tests with y-component only for complex Green tensors

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
def test_vacuum_green_tensor(
    pi_module: PairinteractionModule, use_real: bool, distance_vector_mum: list[float]
) -> None:
    if use_real and distance_vector_mum[1] != 0:
        return  # run tests with y-component only for complex Green tensors

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
    gt: NDArray = (
        np.eye(3) - 3 * np.outer(distance_vector_mum, distance_vector_mum) / distance_mum**2
    ) / distance_au**3
    return gt


def rotation_matrix_from_axis_angle(axis: list[float], angle: float) -> NDArray:
    axis_array = np.array(axis, dtype=float)
    axis_array /= np.linalg.norm(axis_array)
    cross_matrix = np.array(
        [
            [0, -axis_array[2], axis_array[1]],
            [axis_array[2], 0, -axis_array[0]],
            [-axis_array[1], axis_array[0], 0],
        ]
    )
    return np.eye(3) + np.sin(angle) * cross_matrix + (1 - np.cos(angle)) * cross_matrix @ cross_matrix  # type: ignore [no-any-return]


def test_surface_green_tensor_rotation_invariance() -> None:
    rotation = rotation_matrix_from_axis_angle([1, 2, 0.5], np.deg2rad(37))
    pos1 = np.array([0.4, -0.2, 0.8])
    pos2 = np.array([1.1, 0.3, 1.4])
    point_on_plane = np.array([0.0, 0.0, 0.0])

    gt_reference = GreenTensorSurface(
        pos1,
        pos2,
        point_on_plane=point_on_plane,
        surface_normal=[0, 0, 1],
        unit="micrometer",
        static_limit=True,
    ).get(1, 1, transition_energy=0, scaled=True)
    gt_rotated = GreenTensorSurface(
        rotation.T @ pos1,
        rotation.T @ pos2,
        point_on_plane=rotation.T @ point_on_plane,
        surface_normal=rotation.T @ np.array([0, 0, 1]),
        unit="micrometer",
        static_limit=True,
    ).get(1, 1, transition_energy=0, scaled=True)

    np.testing.assert_allclose(
        gt_rotated.m,
        rotation.T @ gt_reference.m @ rotation,
        rtol=1e-10,
        atol=1e-12,
    )


def test_cavity_green_tensor_rotation_invariance() -> None:
    rotation = rotation_matrix_from_axis_angle([0.3, 1.0, -0.2], np.deg2rad(51))
    pos1 = np.array([0.2, -0.1, -0.2])
    pos2 = np.array([0.8, 0.5, 0.4])
    point_on_plane1 = np.array([0.0, 0.0, -1.0])
    point_on_plane2 = np.array([0.0, 0.0, 1.5])

    gt_reference = GreenTensorCavity(
        pos1,
        pos2,
        point_on_plane1=point_on_plane1,
        point_on_plane2=point_on_plane2,
        surface_normal=[0, 0, 1],
        unit="micrometer",
        static_limit=True,
    ).get(1, 1, transition_energy=0, scaled=True)
    gt_rotated = GreenTensorCavity(
        rotation.T @ pos1,
        rotation.T @ pos2,
        point_on_plane1=rotation.T @ point_on_plane1,
        point_on_plane2=rotation.T @ point_on_plane2,
        surface_normal=rotation.T @ np.array([0, 0, 1]),
        unit="micrometer",
        static_limit=True,
    ).get(1, 1, transition_energy=0, scaled=True)

    np.testing.assert_allclose(
        gt_rotated.m,
        rotation.T @ gt_reference.m @ rotation,
        rtol=1e-10,
        atol=1e-12,
    )


def test_surface_green_tensor_rejects_zero_normal() -> None:
    with pytest.raises(ValueError, match="cannot be zero"):
        GreenTensorSurface(
            [0, 0, 0],
            [0, 0, 1],
            point_on_plane=[0, 0, 0],
            surface_normal=[0, 0, 0],
            unit="micrometer",
            static_limit=True,
        )


def test_cavity_green_tensor_rejects_coincident_planes() -> None:
    gt = GreenTensorCavity(
        [0, 0, 0],
        [0, 0, 1],
        point_on_plane1=[0, 0, 0],
        point_on_plane2=[1, 0, 0],
        surface_normal=[0, 0, 1],
        unit="micrometer",
        static_limit=True,
    )

    with pytest.raises(ValueError, match="must be distinct"):
        gt.get(1, 1, transition_energy=0, scaled=True)
