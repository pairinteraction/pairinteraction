# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from collections.abc import Collection, Iterable
from typing import TYPE_CHECKING, Any, Generic, Literal, Optional, TypeVar, Union

import numpy as np
from pint import UnitRegistry
from pint.facets.plain import PlainQuantity
from scipy.sparse import csr_matrix

if TYPE_CHECKING:
    import numpy.typing as npt
    from pint.facets.plain import PlainUnit
    from typing_extensions import Self, TypeAlias

    NDArray: TypeAlias = npt.NDArray[Any]
    PintFloat: TypeAlias = PlainQuantity[float]
    PintArray: TypeAlias = PlainQuantity[NDArray]
    # type ignore here and also below for PlainQuantity[ValueType] because pint has no type support for scipy.csr_matrix
    PintSparse: TypeAlias = PlainQuantity[csr_matrix]  # type: ignore [type-var]
    # and also for complex
    PintComplex: TypeAlias = PlainQuantity[complex]  # type: ignore [type-var]

ureg = UnitRegistry(system="atomic")

Dimension = Literal[
    "electric_field",
    "magnetic_field",
    "distance",
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
    "arbitrary",
    "zero",
]
DimensionLike = Union[Dimension, Iterable[Dimension]]

# some abbreviations: au_time: atomic_unit_of_time; au_current: atomic_unit_of_current; m_e: electron_mass
_CommonUnits: dict[Dimension, str] = {
    "electric_field": "V/cm",  # 1 V/cm = 1.9446903811524456e-10 bohr * m_e / au_current / au_time ** 3
    "magnetic_field": "T",  # 1 T = 4.254382157342044e-06 m_e / au_current / au_time ** 2
    "distance": "micrometer",  # 1 mum = 18897.26124622279 bohr
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
    "arbitrary": "",  # 1 dimensionless
    "zero": "",  # 1 dimensionless
}
BaseUnits: dict[Dimension, "PlainUnit"] = {
    k: ureg.Quantity(1, unit).to_base_units().units for k, unit in _CommonUnits.items()
}
BaseQuantities: dict[Dimension, "PintFloat"] = {k: ureg.Quantity(1, unit) for k, unit in BaseUnits.items()}

Context = Literal["spectroscopy", "Gaussian"]
BaseContexts: dict[Dimension, Context] = {
    "magnetic_field": "Gaussian",
    "energy": "spectroscopy",
}

ValueType = TypeVar("ValueType", bound=Union[float, "NDArray", "csr_matrix"])


class QuantityAbstract(Generic[ValueType]):
    def __init__(self, pint_qty: PlainQuantity[ValueType], dimension: DimensionLike) -> None:  # type: ignore [type-var]
        if not isinstance(pint_qty, ureg.Quantity):
            raise TypeError(f"pint_qty must be a ureg.Quantity, not {type(pint_qty)}")
        self._quantity = pint_qty
        self.dimension: DimensionLike = dimension
        self.check_value_type()

    def check_value_type(self) -> None:
        raise NotImplementedError("This method must be implemented in the derived classes.")

    @classmethod
    def get_base_unit(cls, dimension: DimensionLike) -> str:
        if isinstance(dimension, str):
            return str(BaseUnits[dimension])
        # dimension isinstance Iterable[Dimension]
        return " * ".join(str(BaseUnits[d]) for d in dimension)

    @classmethod
    def get_contexts(cls, dimension: DimensionLike) -> list[Context]:
        if isinstance(dimension, str):
            return [BaseContexts[dimension]] if dimension in BaseContexts else []
        contexts: set[Context] = {BaseContexts[d] for d in dimension if d in BaseContexts}
        return list(contexts)

    @classmethod
    def from_pint(
        cls: "type[Self]",
        value: PlainQuantity[ValueType],  # type: ignore [type-var]
        dimension: DimensionLike,
    ) -> "Self":
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
        cls: "type[Self]",
        value: ValueType,
        unit: str,
        dimension: DimensionLike,
    ) -> "Self":
        """Initialize a Quantity from a value and a unit given as string."""
        if isinstance(value, PlainQuantity):
            raise TypeError("method from_unit: value must be a scalar or an array, not a pint.Quantity")
        return cls(ureg.Quantity(value, unit), dimension)

    @classmethod
    def from_base_unit(
        cls: "type[Self]",
        value: ValueType,
        dimension: DimensionLike,
    ) -> "Self":
        """Initialize a Quantity from a value and a (list of) dimension(s), assume the value is given in base units."""
        unit = cls.get_base_unit(dimension)
        return cls(ureg.Quantity(value, unit), dimension)

    @classmethod
    def from_pint_or_unit(
        cls: "type[Self]",
        value: Union[PlainQuantity[ValueType], ValueType],  # type: ignore [type-var]
        unit: Optional[str],
        dimension: DimensionLike,
    ) -> "Self":
        if unit is None:
            if isinstance(value, PlainQuantity):
                return cls.from_pint(value, dimension)
            if value == 0:
                return cls.from_base_unit(value, dimension)
            raise ValueError("unit must be given if value is not a pint.Quantity")
        assert not isinstance(value, PlainQuantity)
        return cls.from_unit(value, unit, dimension)

    def to_pint(self) -> PlainQuantity[ValueType]:  # type: ignore [type-var]
        """Return the pint.Quantity object."""
        contexts = self.get_contexts(self.dimension)
        base_unit = self.get_base_unit(self.dimension)
        return self._quantity.to(base_unit, *contexts)

    def to_unit(
        self,
        unit: str,
    ) -> ValueType:
        """Return the value of the quantity in the given unit."""
        contexts = self.get_contexts(self.dimension)
        return self._quantity.to(unit, *contexts).magnitude  # type: ignore [no-any-return] # also a problem with pint with sparse matrix

    def to_base_unit(self) -> ValueType:
        """Return the value of the quantity in the base unit."""
        value = self.to_pint().to_base_units()
        return value.magnitude

    def to_pint_or_unit(self, unit: Optional[str]) -> Union[ValueType, PlainQuantity[ValueType]]:  # type: ignore [type-var]
        if unit is None:
            return self.to_pint()
        return self.to_unit(unit)


class QuantityScalar(QuantityAbstract[float]):
    def check_value_type(self) -> None:
        magnitude = self._quantity.magnitude
        if not np.isscalar(magnitude):
            raise TypeError(f"value must be a scalar, not {type(magnitude)}")


class QuantityArray(QuantityAbstract["NDArray"]):
    def check_value_type(self) -> None:
        magnitude = self._quantity.magnitude
        if not isinstance(magnitude, Collection):
            raise TypeError(f"value must be an np.ndarray (or a Collection), not {type(magnitude)}")


class QuantitySparse(QuantityAbstract["csr_matrix"]):
    def check_value_type(self) -> None:
        magnitude = self._quantity.magnitude
        if not isinstance(magnitude, csr_matrix):
            raise TypeError(f"value must be a scipy.sparse.csr_matrix, not {type(magnitude)}")
