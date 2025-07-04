# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later


import numpy as np
import pairinteraction.real as pi
import pytest
from pairinteraction import perturbative
from pairinteraction.units import ureg


@pytest.fixture
def system_pair_sample() -> pi.SystemPair:
    basis = pi.BasisAtom(
        species="Rb",
        n=(59, 63),
        l=(0, 1),
        m=(-1.5, 1.5),
    )
    system = pi.SystemAtom(basis=basis)
    system.set_diamagnetism_enabled(False)
    system.set_magnetic_field([0, 0, 1e-3], "gauss")
    if not system.is_diagonal:
        pi.diagonalize([system], diagonalizer="eigen", sort_by_energy=False)
    basis_pair = pi.BasisPair([system, system])
    system_pair = pi.SystemPair(basis_pair)
    theta = 0
    r = 12
    system_pair.set_distance_vector(r * np.array([np.sin(theta), 0, np.cos(theta)]), "micrometer")
    system_pair.set_interaction_order(3)
    return system_pair


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_c3_with_system(system_pair_sample: pi.SystemPair) -> None:
    """Test whether the C3 coefficient with a given system is calculated correctly."""
    ket_tuple_list = [
        (pi.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5), pi.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5)),
        (pi.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5), pi.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)),
    ]
    c3 = perturbative.get_c3_from_system(
        system_pair=system_pair_sample, ket_tuple_list=ket_tuple_list, unit="planck_constant * gigahertz * micrometer^3"
    )
    assert np.isclose(-0.5 * c3, 3.1515)


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_c3_create_system() -> None:
    """Test whether the C3 coefficient with automatically constructed system is calculated correctly."""
    ket_tuple_list = [
        (pi.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5), pi.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5)),
        (pi.KetAtom("Rb", n=61, l=1, j=1.5, m=0.5), pi.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)),
    ]
    magnetic_field = ureg.Quantity([0, 0, 10], "gauss")
    electric_field = ureg.Quantity([0, 0, 0], "volt/cm")
    distance_vector = ureg.Quantity([0, 0, 500], "micrometer")

    system = perturbative.create_system_for_perturbative(
        ket_tuple_list, electric_field, magnetic_field, distance_vector
    )

    c3 = perturbative.get_c3_from_system(
        system_pair=system, ket_tuple_list=ket_tuple_list, unit="planck_constant * gigahertz * micrometer^3"
    )
    assert np.isclose(-0.5 * c3, 3.2188)


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_c6_with_system(system_pair_sample: pi.SystemPair) -> None:
    """Test whether the C6 coefficient with a given system is calculated correctly."""
    ket_atom = pi.KetAtom(species="Rb", n=61, l=0, j=0.5, m=0.5)
    c6 = perturbative.get_c6_from_system(
        ket_tuple=(ket_atom, ket_atom),
        system_pair=system_pair_sample,
        unit="planck_constant * gigahertz * micrometer^6",
    )
    assert np.isclose(c6, 167.880)


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_c6_create_system() -> None:
    """Test whether the C6 coefficient with automatically constructed system is calculated correctly."""
    magnetic_field = ureg.Quantity([0, 0, 10], "gauss")
    electric_field = ureg.Quantity([0, 0, 0], "volt/cm")
    distance_vector = ureg.Quantity([0, 0, 500], "micrometer")
    ket_atom = pi.KetAtom(species="Rb", n=61, l=0, j=0.5, m=0.5)

    system = perturbative.create_system_for_perturbative(
        [(ket_atom, ket_atom)], electric_field, magnetic_field, distance_vector
    )

    c6 = perturbative.get_c6_from_system(
        ket_tuple=(ket_atom, ket_atom), system_pair=system, unit="planck_constant * gigahertz * micrometer^6"
    )
    assert np.isclose(c6, 169.149)
