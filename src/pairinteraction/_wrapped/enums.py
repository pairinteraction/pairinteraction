from typing import Literal

from pairinteraction import _backend

FloatType = Literal["float32", "float64"]
_FloatTypeDict: dict[FloatType, _backend.FloatType] = {
    "float32": _backend.FloatType.FLOAT32,
    "float64": _backend.FloatType.FLOAT64,
}

OperatorType = Literal[
    "zero",
    "energy",
    "electric_dipole",
    "electric_quadrupole",
    "electric_quadrupole_zero",
    "electric_octupole",
    "magnetic_dipole",
    "arbitrary",
]
_OperatorTypeDict: dict[OperatorType, _backend.OperatorType] = {
    "zero": _backend.OperatorType.ZERO,
    "energy": _backend.OperatorType.ENERGY,
    "electric_dipole": _backend.OperatorType.ELECTRIC_DIPOLE,
    "electric_quadrupole": _backend.OperatorType.ELECTRIC_QUADRUPOLE,
    "electric_quadrupole_zero": _backend.OperatorType.ELECTRIC_QUADRUPOLE_ZERO,
    "electric_octupole": _backend.OperatorType.ELECTRIC_OCTUPOLE,
    "magnetic_dipole": _backend.OperatorType.MAGNETIC_DIPOLE,
    "arbitrary": _backend.OperatorType.ARBITRARY,
}

Parity = Literal["even", "odd", "unknown"]
_ParityDict: dict[Parity, _backend.Parity] = {
    "even": _backend.Parity.EVEN,
    "odd": _backend.Parity.ODD,
    "unknown": _backend.Parity.UNKNOWN,
}


def get_cpp_float_type(float_type: FloatType) -> _backend.FloatType:
    """Convert a python FloatType string to a cpp FloatType enum."""
    if float_type not in _FloatTypeDict:
        raise ValueError(f"Unknown float_type '{float_type}', should be one of {list(_FloatTypeDict.keys())}")
    return _FloatTypeDict[float_type]


def get_cpp_operator_type(operator_type: OperatorType) -> _backend.OperatorType:
    """Convert a python OperatorType string to a cpp OperatorType enum."""
    if operator_type not in _OperatorTypeDict:
        raise ValueError(f"Unknown operator_type '{operator_type}', should be one of {list(_OperatorTypeDict.keys())}")
    return _OperatorTypeDict[operator_type]


def get_cpp_parity(parity: Parity) -> _backend.Parity:
    """Convert a python Parity string to a cpp Parity enum."""
    if parity not in _ParityDict:
        raise ValueError(f"Unknown parity '{parity}', should be one of {list(_ParityDict.keys())}")
    return _ParityDict[parity]
