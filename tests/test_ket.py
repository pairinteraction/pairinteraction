# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import pytest

import pairinteraction.real as pi
from pairinteraction.units import ureg


def test_ket() -> None:
    ket = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    assert ket.species == "Rb"
    assert ket.n == 60
    assert ket.l == 0
    assert pytest.approx(ket.s) == 0.5  # NOSONAR
    assert pytest.approx(ket.j) == 0.5  # NOSONAR
    assert pytest.approx(ket.m) == 0.5  # NOSONAR
    assert ket.parity == "even"
    assert ket.get_energy().units == ureg.Unit("bohr^2 electron_mass/atomic_unit_of_time^2")
    assert pytest.approx(ket.get_energy().magnitude) == 0.15335264334573842  # NOSONAR
    assert pytest.approx(ket.get_energy("GHz")) == 1009011.9215883961  # NOSONAR

    assert ket == pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)

    ket_odd = pi.KetAtom("Rb", n=60, l=1, j=1.5, m=-0.5)
    assert ket_odd.l == 1
    assert pytest.approx(ket_odd.f) == 1.5  # NOSONAR
    assert pytest.approx(ket_odd.m) == -0.5  # NOSONAR
    assert ket_odd.parity == "odd"

    for fmt in ["raw", "ket", "bra"]:
        label = ket_odd.get_label(fmt)
        assert all(str(qn) in label for qn in ["Rb", 60, "P", "3/2", "-1/2"])

    assert ket_odd != ket
    assert ket != pi.KetAtom("Rb", n=60, l=1, j=1.5, m=0.5)

    assert ket.get_matrix_element(ket_odd, "electric_dipole", q=-1) != 0
    assert ket.get_matrix_element(ket_odd, "electric_dipole", q=0) == 0
    assert ket.get_matrix_element(ket_odd, "electric_dipole", q=+1) == 0
