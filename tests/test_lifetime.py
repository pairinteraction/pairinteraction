"""Test calculating lifetimes."""

import numpy as np
import pytest
from scipy.optimize import curve_fit

import pairinteraction.real as pi
from pairinteraction import ureg


def test_lifetime() -> None:
    """Test calculating the lifetime of a state.

    Note that obtaining a reasonable value for the lifetime requires the full database and is not possible with the
    reduced database that is included in the repository!
    """
    ket = pi.KetAtom("Yb174_mqdt", n=60, l=0, j=0, m=0)

    lifetime1 = ket.get_lifetime(temperature=300, temperature_unit="K", unit="us")
    lifetime2 = ket.get_lifetime(temperature=300 * ureg.K, unit="us")
    lifetime3 = ket.get_lifetime(temperature=300 * ureg.K)
    lifetime4 = ket.get_lifetime(unit="us")
    assert lifetime1 == lifetime2 == lifetime3.to(ureg.us).magnitude < lifetime4


@pytest.mark.skip(reason="This test would require the full database.")
def test_lifetime_scaling() -> None:
    """Test the scaling of the lifetime with the principal quantum number."""

    def fit_function(x, a, b):
        return a * x + b

    n_list = list(range(30, 90, 10))

    # S states
    kets = [pi.KetAtom("Rb", n=n, l=0, j=0.5, m=0.5) for n in n_list]
    nu = [ket.nu for ket in kets]
    lifetimes = [ket.get_lifetime(unit="us") for ket in kets]
    popt, _ = curve_fit(fit_function, np.log(nu), np.log(lifetimes))
    assert np.isclose(popt[0], 3, atol=0.1)

    # circular states
    kets = [pi.KetAtom("Rb", n=n, l=n - 1, j=n - 0.5, m=n - 0.5) for n in n_list]
    nu = [ket.nu for ket in kets]
    lifetimes = [ket.get_lifetime(unit="us") for ket in kets]
    popt, _ = curve_fit(fit_function, np.log(nu), np.log(lifetimes))
    assert np.isclose(popt[0], 5, atol=0.1)


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
