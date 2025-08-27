# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from collections.abc import Collection, Iterable
from typing import TYPE_CHECKING, Any, Generic, Literal, Optional, TypeVar, Union

import pint
from pint import UnitRegistry
from pint.facets.plain import PlainQuantity, PlainUnit
from typing_extensions import TypeAlias

if TYPE_CHECKING:
    import numpy.typing as npt
    from pint.facets.numpy.quantity import NumpyQuantity
    from scipy.sparse import csr_matrix

    NDArray: TypeAlias = npt.NDArray[Any]
    ArrayLike: TypeAlias = Union[npt.NDArray[Any], Collection[float]]
    PintFloat: TypeAlias = PlainQuantity[float]
    PintArray: TypeAlias = NumpyQuantity[NDArray]
    PintArrayLike: TypeAlias = Union["PintArray", Collection[Union[float, PintFloat]]]
    # type ignore here and also below for PlainQuantity[ValueType] because pint has no type support for scipy.csr_matrix
    PintSparse: TypeAlias = NumpyQuantity[csr_matrix]  # type: ignore [type-var]
    # and also for complex
    PintComplex: TypeAlias = PlainQuantity[complex]  # type: ignore [type-var]

ureg = UnitRegistry(system="atomic")

Dimension = Literal[
    "electric_field",
    "magnetic_field",
    "distance",
    "inverse_distance",
    "energy",
    "charge",
    "velocity",
    "temperature",
    "time",
    "transition_rate",
    "electric_dipole",
    "electric_quadrupole",
    "electric_quadrupole_zero",
    "electric_octupole",
    "magnetic_dipole",
    "c3",
    "c6",
    "green_tensor_00",
    "arbitrary",
    "zero",
]
DimensionLike = Union[Dimension, Iterable[Dimension]]

# some abbreviations: au_time: atomic_unit_of_time; au_current: atomic_unit_of_current; m_e: electron_mass
_CommonUnits: dict[Dimension, str] = {
    "electric_field": "V/cm",  # 1 V/cm = 1.9446903811524456e-10 bohr * m_e / au_current / au_time ** 3
    "magnetic_field": "T",  # 1 T = 4.254382157342044e-06 m_e / au_current / au_time ** 2
    "distance": "micrometer",  # 1 mum = 18897.26124622279 bohr
    "inverse_distance": "1 / micrometer",
    "energy": "hartree",  # 1 hartree = 1 bohr ** 2 * m_e / au_time ** 2
    "charge": "e",  # 1 e = 1 au_current * au_time
    "velocity": "speed_of_light",  # 1 c = 137.03599908356244 bohr / au_time
    "temperature": "K",  # 1 K = 3.1668115634555572e-06 atomic_unit_of_temperature
    "time": "s",  # 1 s = 4.134137333518244e+16 au_time
    "transition_rate": "1/s",  # 1 / s = 2.4188843265856806e-17 * 1 / au_time
    "electric_dipole": "e * a0",  # 1 e * a0 = 1 au_current * au_time * bohr
    "electric_quadrupole": "e * a0^2",  # 1 e * a0^2 = 1 au_current * au_time * bohr ** 2
    "electric_quadrupole_zero": "e * a0^2",  # 1 e * a0^2 = 1 au_current * au_time * bohr ** 2
    "electric_octupole": "e * a0^3",  # 1 e * a0^3 = 1 au_current * au_time * bohr ** 3
    "magnetic_dipole": "hbar e / m_e",  # 1 hbar e / m_e = 1 au_current * bohr ** 2
    "c3": "hartree * bohr^3",  # 1 hartree * bohr^3 = 1 bohr ** 3 * m_e / au_time ** 2
    "c6": "hartree * bohr^6",  # 1 hartree * bohr^6 = 1 bohr ** 6 * m_e / au_time ** 2
    "green_tensor_00": "hartree / e^2",  # unit for green tensor with kappa1 = kappa2 = 0
    "arbitrary": "",  # 1 dimensionless
    "zero": "",  # 1 dimensionless
}
AtomicUnits: dict[Dimension, "PlainUnit"] = {
    k: ureg.Quantity(1, unit).to_base_units().units for k, unit in _CommonUnits.items()
}

Context = Literal["spectroscopy", "Gaussian"]
BaseContexts: dict[Dimension, Context] = {
    "magnetic_field": "Gaussian",
    "energy": "spectroscopy",
    "c3": "spectroscopy",
    "c6": "spectroscopy",
}

ValueType = TypeVar("ValueType", bound=Union[float, "NDArray", "csr_matrix"])
PintType = TypeVar("PintType", bound=Union["PintFloat", "PintArray", "PintSparse"])  # type: ignore [type-var]
ValueTypeLike = TypeVar("ValueTypeLike", bound=Union[float, "ArrayLike", "csr_matrix"])
PintTypeLike = TypeVar("PintTypeLike", bound=Union["PintFloat", "PintArrayLike", "PintSparse"])
T = TypeVar("T", bound=PlainQuantity[Any])


class UnitConverterGeneric(Generic[ValueType, PintType, ValueTypeLike, PintTypeLike]):
    @staticmethod
    def _get_unit_au_from_dimension(dimension: DimensionLike) -> str:
        """Get the unit for a given dimension."""
        if isinstance(dimension, str):
            return str(AtomicUnits[dimension])
        return " * ".join(str(AtomicUnits[d]) for d in dimension)

    @staticmethod
    def _get_contexts_from_dimension(dimension: DimensionLike) -> list[Context]:
        if isinstance(dimension, str):
            return [BaseContexts[dimension]] if dimension in BaseContexts else []
        contexts: set[Context] = {BaseContexts[d] for d in dimension if d in BaseContexts}
        return list(contexts)

    @staticmethod
    def _pint_to_pint(
        quantity: T,
        unit: str,
        dimension: DimensionLike,
    ) -> T:
        """Convert a pint.Quantity to another pint.Quantity with a different unit."""
        contexts = UnitConverterGeneric._get_contexts_from_dimension(dimension)
        try:
            return quantity.to(unit, *contexts)  # type: ignore [return-value]
        except pint.errors.DimensionalityError:
            # pint uses e.g. the context "spectroscopy" to convert "hartree" -> "GHz"
            # however, something like "hartree * bohr^3" -> "GHz * bohr^3" does not work
            # the following is a workaround for this kind of conversions
            if "spectroscopy" in contexts:
                quantity = quantity * ureg.Quantity(1, "GHz") / ureg.Quantity(1, "GHz").to("hartree", "spectroscopy")
                return quantity.to(unit, *contexts)  # type: ignore [return-value]
            raise

    @classmethod
    def au_to_pint(
        cls,
        value: ValueTypeLike,
        dimension: DimensionLike,
    ) -> PintType:
        """Convert a value in atomic units to a pint.Quantity."""
        unit_au = cls._get_unit_au_from_dimension(dimension)
        return cls.user_to_pint(value, unit_au, dimension)

    @classmethod
    def user_to_pint(
        cls,
        qty_val: Union[PintTypeLike, ValueTypeLike],
        unit: Optional[str],
        dimension: DimensionLike,
    ) -> PintType:
        """Convert a user-defined value to a pint.Quantity."""
        unit_au = cls._get_unit_au_from_dimension(dimension)
        quantity: PintType
        if isinstance(qty_val, PlainQuantity):
            if unit is not None:
                raise ValueError("unit must be None if qty_val is a pint.Quantity")
            quantity = qty_val  # type: ignore [assignment]
        else:
            if isinstance(qty_val, Iterable) and unit is None:
                unit = unit_au
                new_qty_val = []
                for qv in qty_val:
                    if qv == 0:
                        new_qty_val.append(0)
                    elif isinstance(qv, PlainQuantity):
                        value = cls._pint_to_pint(qv, unit_au, dimension).magnitude
                        new_qty_val.append(value)
                    else:
                        raise ValueError("unit must be given if qty_val is not a pint.Quantity and non-zero")
                qty_val = new_qty_val  # type: ignore [assignment]
            elif unit is None:
                unit = unit_au
                if qty_val != 0:
                    raise ValueError("unit must be given if qty_val is not a pint.Quantity")
            quantity = ureg.Quantity(qty_val, unit)  # type: ignore [assignment]
        return cls._pint_to_pint(quantity, unit_au, dimension)

    @classmethod
    def pint_to_au(
        cls,
        quantity: PintType,
        dimension: DimensionLike,
    ) -> ValueType:
        """Convert a pint.Quantity to a value in atomic units."""
        unit_au = cls._get_unit_au_from_dimension(dimension)
        return cls._pint_to_pint(quantity, unit_au, dimension).magnitude  # type: ignore [return-value]

    @classmethod
    def user_to_au(
        cls,
        qty_val: Union[PintTypeLike, ValueTypeLike],
        unit: Optional[str],
        dimension: DimensionLike,
    ) -> ValueType:
        """Convert a user-defined value to a value in atomic units."""
        quantity = cls.user_to_pint(qty_val, unit, dimension)
        return cls.pint_to_au(quantity, dimension)

    @classmethod
    def au_to_user(
        cls,
        value: ValueTypeLike,
        dimension: DimensionLike,
        to_unit: Optional[str],
    ) -> Union[PintType, ValueType]:
        """Convert a value in atomic units to a user-defined value."""
        quantity = cls.au_to_pint(value, dimension)
        return cls.pint_to_user(quantity, dimension, to_unit)

    @classmethod
    def pint_to_user(
        cls,
        quantity: PintType,
        dimension: DimensionLike,
        to_unit: Optional[str],
    ) -> Union[PintType, ValueType]:
        """Convert a pint.Quantity to a user-defined value."""
        if to_unit is None:
            unit_au = cls._get_unit_au_from_dimension(dimension)
            return cls._pint_to_pint(quantity, unit_au, dimension)
        return cls._pint_to_pint(quantity, to_unit, dimension).magnitude  # type: ignore [return-value]


class UnitConverterScalar(UnitConverterGeneric[float, "PintFloat", float, "PintFloat"]):
    """Unit converter for scalar values."""


class UnitConverterArray(UnitConverterGeneric["NDArray", "PintArray", "ArrayLike", "PintArrayLike"]):
    """Unit converter for array values."""


class UnitConverterSparse(UnitConverterGeneric["csr_matrix", "PintSparse", "csr_matrix", "PintSparse"]):
    """Unit converter for sparse matrix values."""
