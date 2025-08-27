# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import Literal

import pairinteraction.real as pi
import pytest
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

    formats: list[Literal["raw", "ket", "bra"]] = ["raw", "ket", "bra"]
    for fmt in formats:
        label = ket_odd.get_label(fmt)
        assert all(str(qn) in label for qn in ["Rb", 60, "P", "3/2", "-1/2"])

    assert ket_odd != ket
    assert ket != pi.KetAtom("Rb", n=60, l=1, j=1.5, m=0.5)

    assert ket.get_matrix_element(ket_odd, "electric_dipole", q=-1) != 0
    assert ket.get_matrix_element(ket_odd, "electric_dipole", q=0) == 0
    assert ket.get_matrix_element(ket_odd, "electric_dipole", q=+1) == 0


def test_ket_equal() -> None:
    ket1 = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket2 = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket3 = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=-0.5)
    ket4 = pi.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    assert ket1 == ket2
    assert ket1 != ket3
    assert ket1 != ket4

    ket1 = pi.KetAtom("Sr88_singlet", n=60, l=1, j=1, m=0)
    ket2 = pi.KetAtom("Sr88_triplet", n=60, l=1, j=1, m=0)
    assert ket1 != ket2
