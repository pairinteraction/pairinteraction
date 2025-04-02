# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Tests for atomic unit conversions using pint."""

import pytest
from pint import UnitRegistry

from tests.constants import (
    GAUSS_IN_ATOMIC_UNITS,
    HARTREE_IN_GHZ,
    HARTREE_IN_INVERSE_CM,
    HARTREE_IN_JOULES,
    VOLT_PER_CM_IN_ATOMIC_UNITS,
)


def test_hartree_to_joules(ureg: UnitRegistry) -> None:
    """Test conversion from Hartree to Joules."""
    one_hartree = 1 * ureg.hartree
    one_hartree_in_joules = one_hartree.to("joule")
    assert pytest.approx(one_hartree_in_joules.magnitude, rel=1e-12) == HARTREE_IN_JOULES  # NOSONAR


def test_joules_to_hartree(ureg: UnitRegistry) -> None:
    """Test conversion from Joules to Hartree."""
    one_hartree_in_joules = HARTREE_IN_JOULES * ureg.joule
    result = one_hartree_in_joules.to_base_units()
    assert pytest.approx(result.magnitude, rel=1e-12) == 1.0  # NOSONAR


def test_hartree_to_thz(ureg: UnitRegistry) -> None:
    """Test conversion from Hartree to THz."""
    one_hartree = 1 * ureg.hartree
    one_hartree_in_ghz = one_hartree.to("gigahertz", "spectroscopy")
    assert pytest.approx(one_hartree_in_ghz.magnitude, rel=1e-12) == HARTREE_IN_GHZ  # NOSONAR


def test_hartree_to_inverse_cm(ureg: UnitRegistry) -> None:
    """Test conversion from Hartree to inverse cm."""
    one_hartree = 1 * ureg.hartree
    one_hartree_in_inverse_cm = one_hartree.to("1/cm", "spectroscopy")
    assert pytest.approx(one_hartree_in_inverse_cm.magnitude, rel=1e-12) == HARTREE_IN_INVERSE_CM  # NOSONAR


def test_electric_field_to_atomic_units(ureg: UnitRegistry) -> None:
    """Test conversion from V/cm to atomic units of electric field."""
    one_v_per_cm = 1 * ureg.volt / ureg.centimeter
    one_v_per_cm_in_atomic_units = one_v_per_cm.to_base_units()
    assert pytest.approx(one_v_per_cm_in_atomic_units.magnitude, rel=1e-12) == VOLT_PER_CM_IN_ATOMIC_UNITS  # NOSONAR


def test_magnetic_field_to_atomic_units(ureg: UnitRegistry) -> None:
    """Test conversion from Gauss to atomic units of magnetic field."""
    one_gauss = 1e-4 * ureg.tesla
    one_gauss_in_atomic_units = one_gauss.to_base_units()
    assert pytest.approx(one_gauss_in_atomic_units.magnitude, rel=1e-12) == GAUSS_IN_ATOMIC_UNITS  # NOSONAR
