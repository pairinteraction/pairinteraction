"""
Tests for atomic unit conversions using pint.
"""

import pytest
from pint import UnitRegistry

HARTREE_IN_JOULES = 4.3597447222060e-18
HARTREE_IN_THZ = 6579.683920502
HARTREE_IN_INVERSE_CM = 219474.63136320


def test_hartree_to_joules(ureg: UnitRegistry) -> None:
    """Test conversion from Hartree to Joules."""
    one_hartree = 1 * ureg.hartree
    one_hartree_in_joules = one_hartree.to("joule")
    assert pytest.approx(one_hartree_in_joules.magnitude, rel=1e-12) == HARTREE_IN_JOULES


def test_joules_to_hartree(ureg: UnitRegistry) -> None:
    """Test conversion from Joules to Hartree."""
    one_hartree_in_joules = HARTREE_IN_JOULES * ureg.joule
    result = one_hartree_in_joules.to_base_units()
    assert pytest.approx(result.magnitude, rel=1e-12) == 1.0


def test_hartree_to_thz(ureg: UnitRegistry) -> None:
    """Test conversion from Hartree to THz."""
    one_hartree = 1 * ureg.hartree
    one_hartree_in_thz = one_hartree.to("terahertz", "spectroscopy")
    assert pytest.approx(one_hartree_in_thz.magnitude, rel=1e-12) == HARTREE_IN_THZ


def test_hartree_to_inverse_cm(ureg: UnitRegistry) -> None:
    """Test conversion from Hartree to inverse cm."""
    one_hartree = 1 * ureg.hartree
    one_hartree_in_inverse_cm = one_hartree.to("1/cm", "spectroscopy")
    assert pytest.approx(one_hartree_in_inverse_cm.magnitude, rel=1e-12) == HARTREE_IN_INVERSE_CM


def test_electric_field_to_atomic_units(ureg: UnitRegistry) -> None:
    """Test conversion from V/cm to atomic units of electric field."""
    one_v_per_cm = 1 * ureg.volt / ureg.centimeter
    one_v_per_cm_in_atomic_units = one_v_per_cm.to_base_units()
    assert pytest.approx(one_v_per_cm_in_atomic_units.magnitude, rel=1e-12) == 1 / 5.14220675112e9


def test_magnetic_field_to_atomic_units(ureg: UnitRegistry) -> None:
    """Test conversion from Gauss to atomic units of magnetic field."""
    one_gauss = 1e-4 * ureg.tesla
    one_gauss_in_atomic_units = one_gauss.to_base_units()
    assert pytest.approx(one_gauss_in_atomic_units.magnitude, rel=1e-12) == 1 / 2.35051757077e9
