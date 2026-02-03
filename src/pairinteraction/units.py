# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from collections.abc import Collection, Iterable
from typing import TYPE_CHECKING, Any, Generic, Literal, TypeVar, Union

import numpy as np
import pint
from pint import UnitRegistry
from pint.facets.plain import PlainQuantity
from scipy.sparse import csr_matrix

if TYPE_CHECKING:
    import numpy.typing as npt
    from pint.facets.plain import PlainUnit
    from typing_extensions import Self, TypeAlias

    NDArray: TypeAlias = npt.NDArray[Any]
    ArrayLike: TypeAlias = Union[npt.NDArray[Any], Collection[float]]
    PintFloat: TypeAlias = PlainQuantity[float]
    PintArray: TypeAlias = PlainQuantity[NDArray]
    PintArrayLike: TypeAlias = Union["PintArray", Collection[Union[float, "PintFloat"]]]
    # type ignore here and also below for PlainQuantity[ValueType] because pint has no type support for scipy.csr_matrix
    PintSparse: TypeAlias = PlainQuantity[csr_matrix]  # type: ignore [type-var]
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
    "green_tensor_dd",
    "scaled_green_tensor_dd",
    "identity",
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
    "green_tensor_00": "meter",  # unit for green tensor with kappa1 = kappa2 = 0
    "green_tensor_dd": "1 / meter",  # unit for green tensor with kappa1 = kappa2 = 1
    "scaled_green_tensor_dd": "hartree / (e * meter)^2",  # unit for scaled green tensor with kappa1 = kappa2 = 1
    "identity": "",  # 1 dimensionless
    "arbitrary": "",  # 1 dimensionless
    "zero": "",  # 1 dimensionless
}
AtomicUnits: dict[Dimension, PlainUnit] = {
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
ValueTypeLike = TypeVar("ValueTypeLike", bound=Union[float, "ArrayLike", "csr_matrix"])


class QuantityAbstract(Generic[ValueTypeLike, ValueType]):
    def __init__(self, pint_qty: PlainQuantity[ValueType], dimension: DimensionLike) -> None:  # type: ignore [type-var]
        if not isinstance(pint_qty, ureg.Quantity):
            raise TypeError(f"pint_qty must be a ureg.Quantity, not {type(pint_qty)}")
        self._quantity = pint_qty
        self.dimension: DimensionLike = dimension
        self.check_value_type()

    def check_value_type(self) -> None:
        raise NotImplementedError("This method must be implemented in the derived classes.")

    @classmethod
    def get_atomic_unit(cls, dimension: DimensionLike) -> str:
        if isinstance(dimension, str):
            return str(AtomicUnits[dimension])
        # dimension isinstance Iterable[Dimension]
        return " * ".join(str(AtomicUnits[d]) for d in dimension)

    @classmethod
    def get_contexts(cls, dimension: DimensionLike) -> list[Context]:
        if isinstance(dimension, str):
            return [BaseContexts[dimension]] if dimension in BaseContexts else []
        contexts: set[Context] = {BaseContexts[d] for d in dimension if d in BaseContexts}
        return list(contexts)

    @classmethod
    def from_pint(
        cls: type[Self],
        value: PlainQuantity[ValueType],  # type: ignore [type-var]
        dimension: DimensionLike,
    ) -> Self:
        """Initialize a Quantity from a ureg.Quantity."""
        if isinstance(value, ureg.Quantity):
            return cls(value, dimension)
        if isinstance(value, PlainQuantity):
            raise TypeError(
                "Only use pint quantities genereated by pairinteraction.ureg and not by a different pint.UnitRegistry."
            )
        raise ValueError("method from_pint: value must be a pint.Quantity")

    @classmethod
    def from_unit(
        cls: type[Self],
        value: ValueTypeLike,
        unit: str,
        dimension: DimensionLike,
    ) -> Self:
        """Initialize a Quantity from a value and a unit given as string."""
        if isinstance(value, PlainQuantity):
            raise TypeError("method from_unit: value must be a scalar or an array, not a pint.Quantity")
        return cls(ureg.Quantity(value, unit), dimension)

    @classmethod
    def from_au(
        cls: type[Self],
        value: ValueTypeLike,
        dimension: DimensionLike,
    ) -> Self:
        """Initialize a Quantity from a value in atomic units (a.u.) and a (list of) dimension(s)."""
        unit = cls.get_atomic_unit(dimension)
        return cls(ureg.Quantity(value, unit), dimension)

    @classmethod
    def from_pint_or_unit(
        cls: type[Self],
        value: PlainQuantity[ValueType] | ValueTypeLike,  # type: ignore [type-var]
        unit: str | None,
        dimension: DimensionLike,
    ) -> Self:
        if unit is None:
            if isinstance(value, PlainQuantity):
                return cls.from_pint(value, dimension)
            if np.all(value == 0):
                return cls.from_au(value, dimension)
            raise ValueError("unit must be given if value is not a pint.Quantity")
        assert not isinstance(value, PlainQuantity)
        return cls.from_unit(value, unit, dimension)

    def to_pint(self) -> PlainQuantity[ValueType]:  # type: ignore [type-var]
        """Return the pint.Quantity object."""
        contexts = self.get_contexts(self.dimension)
        atomic_unit = self.get_atomic_unit(self.dimension)
        return self._quantity.to(atomic_unit, *contexts)

    def to_unit(
        self,
        unit: str,
    ) -> ValueType:
        """Return the value of the quantity in the given unit."""
        contexts = self.get_contexts(self.dimension)
        try:
            return self._quantity.to(unit, *contexts).magnitude  # type: ignore [no-any-return] # also a problem with pint with sparse matrix
        except pint.errors.DimensionalityError:
            # pint uses e.g. the context "spectroscopy" to convert "hartree" -> "GHz"
            # however, something like "hartree * bohr^3" -> "GHz * bohr^3" does not work
            # the following is a workaround for this kind of conversions
            if "spectroscopy" in contexts:
                q = self._quantity * ureg.Quantity(1, "GHz") / ureg.Quantity(1, "GHz").to("hartree", "spectroscopy")
                return q.to(unit, *contexts).magnitude  # type: ignore [no-any-return]
            raise

    def to_au(self) -> ValueType:
        """Return the value of the quantity in atomic units (a.u.)."""
        value = self.to_pint().to_base_units()
        return value.magnitude

    def to_pint_or_unit(self, unit: str | None) -> ValueType | PlainQuantity[ValueType]:  # type: ignore [type-var]
        if unit is None:
            return self.to_pint()
        return self.to_unit(unit)

    @classmethod
    def convert_user_to_au(
        cls,
        value: PlainQuantity[ValueType] | ValueTypeLike,  # type: ignore [type-var]
        unit: str | None,
        dimension: DimensionLike,
    ) -> ValueType:
        return cls.from_pint_or_unit(value, unit, dimension).to_au()

    @classmethod
    def convert_user_to_pint(
        cls,
        value: PlainQuantity[ValueType] | ValueTypeLike,  # type: ignore [type-var]
        unit: str | None,
        dimension: DimensionLike,
    ) -> PlainQuantity[ValueType]:  # type: ignore [type-var]
        return cls.from_pint_or_unit(value, unit, dimension).to_pint()

    @classmethod
    def convert_au_to_user(
        cls,
        values_au: ValueTypeLike,
        dimension: DimensionLike,
        unit: str | None,
    ) -> ValueType | PlainQuantity[ValueType]:  # type: ignore [type-var]
        return cls.from_au(values_au, dimension).to_pint_or_unit(unit)

    @classmethod
    def convert_pint_to_user(
        cls,
        value_pint: PlainQuantity[ValueType],  # type: ignore [type-var]
        dimension: DimensionLike,
        unit: str | None,
    ) -> ValueType | PlainQuantity[ValueType]:  # type: ignore [type-var]
        return cls.from_pint(value_pint, dimension).to_pint_or_unit(unit)


class QuantityScalar(QuantityAbstract[float, float]):
    def check_value_type(self) -> None:
        magnitude = self._quantity.magnitude
        if not np.isscalar(magnitude):
            raise TypeError(f"value must be a scalar, not {type(magnitude)}")


class QuantityArray(QuantityAbstract["ArrayLike", "NDArray"]):
    def check_value_type(self) -> None:
        magnitude = self._quantity.magnitude
        if not isinstance(magnitude, Collection):
            raise TypeError(f"value must be an np.ndarray (or a Collection), not {type(magnitude)}")


class QuantitySparse(QuantityAbstract["csr_matrix", "csr_matrix"]):
    def check_value_type(self) -> None:
        magnitude = self._quantity.magnitude
        if not isinstance(magnitude, csr_matrix):
            raise TypeError(f"value must be a scipy.sparse.csr_matrix, not {type(magnitude)}")
