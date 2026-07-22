# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import pytest
from pairinteraction.units import ureg

if TYPE_CHECKING:
    from .utils import PairinteractionModule


def test_ket(pi_module: PairinteractionModule) -> None:
    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
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

    assert ket == pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)

    ket_odd = pi_module.KetAtom("Rb", n=60, l=1, j=1.5, m=-0.5)
    assert ket_odd.l == 1
    assert pytest.approx(ket_odd.f) == 1.5  # NOSONAR
    assert pytest.approx(ket_odd.m) == -0.5  # NOSONAR
    assert ket_odd.parity == "odd"

    formats: list[Literal["raw", "ket", "bra"]] = ["raw", "ket", "bra"]
    for fmt in formats:
        label = ket_odd.get_label(fmt)
        assert all(str(qn) in label for qn in ["Rb", 60, "P", "3/2", "-1/2"])

    assert ket_odd != ket
    assert ket != pi_module.KetAtom("Rb", n=60, l=1, j=1.5, m=0.5)

    assert ket.get_matrix_element(ket_odd, "electric_dipole", q=-1) != 0
    assert ket.get_matrix_element(ket_odd, "electric_dipole", q=0) == 0
    assert ket.get_matrix_element(ket_odd, "electric_dipole", q=+1) == 0


def test_get_label_mqdt(pi_module: PairinteractionModule) -> None:
    ket1 = pi_module.KetAtom("Yb171_mqdt", nu=55.5, l=0, f=1.5, m=1.5)
    assert ket1.get_label("raw") == "Yb171:S=1.0,nu=55.6,L=0.0,F=3/2,3/2"
    ket2 = pi_module.KetAtom("Yb171_mqdt", nu=55.1, l=1, f=2.5, m=2.5)
    assert ket2.get_label("raw") == "Yb171:S=1.0,nu=55.1,L=1.0,F=5/2,5/2"
    ket3 = pi_module.KetAtom("Yb174_mqdt", nu=60, l=1, f=1, m=1)
    assert ket3.get_label("raw") == "Yb174:S=0.0,nu=60.0,L=1.0,J=1,1"


def test_ket_equal(pi_module: PairinteractionModule) -> None:
    ket1 = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket2 = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket3 = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=-0.5)
    ket4 = pi_module.KetAtom("Rb", n=61, l=0, j=0.5, m=0.5)
    assert ket1 == ket2
    assert ket1 != ket3
    assert ket1 != ket4

    ket1 = pi_module.KetAtom("Sr88_singlet", n=60, l=1, j=1, m=0)
    ket2 = pi_module.KetAtom("Sr88_triplet", n=60, l=1, j=1, m=0)
    assert ket1 != ket2


def test_ket_to_state(pi_module: PairinteractionModule) -> None:
    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    state = ket.to_state()

    # the state is trivial, i.e. a single ket with a coefficient of one
    assert isinstance(state, pi_module.StateAtom)
    assert state.number_of_kets == 1
    assert state.is_canonical
    assert state.is_normalized()
    assert state.get_corresponding_ket() == ket
    assert pytest.approx(state.get_coefficients()) == [1.0]  # NOSONAR
    assert pytest.approx(state.get_overlap(ket)) == 1.0  # NOSONAR

    # the state has the same data type (real/complex) as the ket, so both can be combined
    combined = (state + pi_module.KetAtom("Rb", n=60, l=1, j=1.5, m=0.5).to_state()).normalize()
    assert combined.number_of_kets == 2
    assert pytest.approx(combined.get_overlap(ket)) == 0.5  # NOSONAR
