# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Test calculating lifetimes."""

from typing import TYPE_CHECKING

import numpy as np
import pairinteraction.real as pi
import pytest
from pairinteraction import ureg
from scipy.optimize import curve_fit

if TYPE_CHECKING:
    from pairinteraction.units import NDArray


def test_lifetime() -> None:
    """Test calculating the lifetime of a state.

    Note that obtaining a reasonable value for the lifetime requires the full database and is not possible with the
    reduced database that is included in the repository!
    """
    ket = pi.KetAtom("Yb174_mqdt", n=64, l=0, j=0, m=0)

    # We use n=64 here, because this corresponds to nu~60 and we include
    # states with 50<nu<70 in the limited database that comes with the repository. In
    # addition, states with nu<20 are included so that the decay to the ground state
    # is also captured. Note that the relative error of the calculated
    # lifetime versus the actual lifetime is still quite large because of the limited
    # number of states. The full database is used and accurate results can be obtained
    # if the test is run with `pytest --database-dir "" --download-missing`.

    lifetime1 = ket.get_lifetime(temperature=300, temperature_unit="K", unit="us")
    lifetime2 = ket.get_lifetime(temperature=300 * ureg.K, unit="us")
    lifetime3 = ket.get_lifetime(temperature=300 * ureg.K)
    lifetime4 = ket.get_lifetime(unit="us")

    assert lifetime1 == lifetime2 == lifetime3.to(ureg.us).magnitude < lifetime4
    assert pytest.approx(lifetime1, rel=0.15) == 142.04845576112646  # NOSONAR
    assert pytest.approx(lifetime4, rel=0.15) == 494.1653414977515  # NOSONAR


def test_lifetime_scaling() -> None:
    """Test the scaling of the lifetime with the principal quantum number."""

    def fit_function(x: "NDArray", a: float, b: float) -> "NDArray":
        return a * x + b

    n_list = list(range(60, 70, 1))

    # S states
    kets = [pi.KetAtom("Rb", n=n, l=0, j=0.5, m=0.5) for n in n_list]
    nu = [ket.nu for ket in kets]
    lifetimes = [ket.get_lifetime(unit="us") for ket in kets]
    popt, _ = curve_fit(fit_function, np.log(nu), np.log(lifetimes))  # type: ignore [arg-type]
    assert np.isclose(popt[0], 3, atol=0.02)

    # Circular states
    try:
        kets = [pi.KetAtom("Rb", n=n, l=n - 1, j=n - 0.5, m=n - 0.5) for n in n_list]
    except ValueError as err:
        # If the limited database which comes with the repository is used, creating the
        # kets will fail because the circular states are not included in the database.
        # This is expected.
        if "No state found" not in str(err):
            raise
    else:
        nu = [ket.nu for ket in kets]
        lifetimes = [ket.get_lifetime(unit="us") for ket in kets]
        popt, _ = curve_fit(fit_function, np.log(nu), np.log(lifetimes))  # type: ignore [arg-type]
        assert np.isclose(popt[0], 5, atol=0.02)


def test_transition_rates() -> None:
    """Test calculating transition rates to other states.

    Note that obtaining a reasonable value for the transition rates requires the full database and is not possible
    with the reduced database that is included in the repository!
    """
    ket = pi.KetAtom("Yb174_mqdt", n=60, l=0, j=0, m=0)

    lifetime = ket.get_lifetime(temperature=300, temperature_unit="K", unit="us")
    kets_sp, rates_sp = ket.get_spontaneous_transition_rates(unit="MHz")
    kets_bbr, rates_bbr = ket.get_black_body_transition_rates(temperature=300, temperature_unit="K", unit="MHz")

    assert len(rates_sp) == len(kets_sp)
    assert len(rates_bbr) == len(kets_bbr)
    assert np.isclose(1 / (sum(rates_sp) + sum(rates_bbr)), lifetime)
