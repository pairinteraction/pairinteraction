from collections.abc import Collection, Iterable
from typing import TYPE_CHECKING, Generic, Literal, Optional, TypeVar, Union

import numpy as np
from pint import UnitRegistry
from pint.facets.plain import PlainQuantity
from scipy.sparse import csr_matrix

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from pint.facets.plain import PlainUnit
    from typing_extensions import Self

ureg = UnitRegistry(system="atomic")

Dimension = Literal[
    "ELECTRIC_FIELD",
    "MAGNETIC_FIELD",
    "DISTANCE",
    "ENERGY",
    "CHARGE",
    "VELOCITY",
    "TEMPERATURE",
    "TIME",
    "TRANSITION_RATE",
    "ELECTRIC_DIPOLE",
    "ELECTRIC_QUADRUPOLE",
    "ELECTRIC_QUADRUPOLE_ZERO",
    "ELECTRIC_OCTUPOLE",
    "MAGNETIC_DIPOLE",
    "C3",
    "C6",
    "ARBITRARY",
    "ZERO",
]
DimensionLike = Union[Dimension, Iterable[Dimension]]

# au_time = atomic_unit_of_time; au_current = atomic_unit_of_current; m_e = electron_mass
_CommonUnits: dict[Dimension, str] = {
    "ELECTRIC_FIELD": "V/cm",  # 1 V/cm = 1.9446903811524456e-10 bohr * m_e / au_current / au_time ** 3
    "MAGNETIC_FIELD": "T",  # 1 T = 4.254382157342044e-06 m_e / au_current / au_time ** 2
    "DISTANCE": "micrometer",  # 1 mum = 18897.26124622279 bohr
    "ENERGY": "hartree",  # 1 hartree = 1 bohr ** 2 * m_e / au_time ** 2
    "CHARGE": "e",  # 1 e = 1 au_current * au_time
    "VELOCITY": "speed_of_light",  # 1 c = 137.03599908356244 bohr / au_time
    "TEMPERATURE": "K",  # 1 K = 3.1668115634555572e-06 atomic_unit_of_temperature
    "TIME": "s",  # 1 s = 4.134137333518244e+16 au_time
    "TRANSITION_RATE": "1/s",  # 1 / s = 2.4188843265856806e-17 * 1 / au_time
    "ELECTRIC_DIPOLE": "e * a0",  # 1 e * a0 = 1 au_current * au_time * bohr
    "ELECTRIC_QUADRUPOLE": "e * a0^2",  # 1 e * a0^2 = 1 au_current * au_time * bohr ** 2
    "ELECTRIC_QUADRUPOLE_ZERO": "e * a0^2",  # 1 e * a0^2 = 1 au_current * au_time * bohr ** 2
    "ELECTRIC_OCTUPOLE": "e * a0^3",  # 1 e * a0^3 = 1 au_current * au_time * bohr ** 3
    "MAGNETIC_DIPOLE": "hbar e / m_e",  # 1 hbar e / m_e = 1 au_current * bohr ** 2
    "C3": "hartree * bohr^3",  # 1 hartree * bohr^3 = 1 bohr ** 3 * m_e / au_time ** 2
    "C6": "hartree * bohr^6",  # 1 hartree * bohr^6 = 1 bohr ** 6 * m_e / au_time ** 2
    "ARBITRARY": "",  # 1 dimensionless
    "ZERO": "",  # 1 dimensionless
}
BaseUnits: dict[Dimension, "PlainUnit"] = {
    k: ureg.Quantity(1, unit).to_base_units().units for k, unit in _CommonUnits.items()
}
BaseQuantities: dict[Dimension, "PlainQuantity[float]"] = {k: ureg.Quantity(1, unit) for k, unit in BaseUnits.items()}

Context = Literal["spectroscopy", "Gaussian"]
BaseContexts: dict[Dimension, Context] = {
    "MAGNETIC_FIELD": "Gaussian",
    "ENERGY": "spectroscopy",
}

ValueType = TypeVar("ValueType", bound=Union[float, "NDArray", "csr_matrix"])


class QuantityAbstract(Generic[ValueType]):
    def __init__(self, pint_qty: PlainQuantity[ValueType], dimension: DimensionLike) -> None:
        if not isinstance(pint_qty, ureg.Quantity):
            raise ValueError(f"pint_qty must be a ureg.Quantity, not {type(pint_qty)}")
        self._quantity = pint_qty
        self.dimension: DimensionLike = dimension
        self.check_value_type()

    def check_value_type(self) -> None:
        raise NotImplementedError("This method must be implemented in the derived classes.")

    @classmethod
    def get_base_unit(cls, dimension: DimensionLike) -> str:
        if isinstance(dimension, str):
            return str(BaseUnits[dimension])
        else:  # dimension isinstance Iterable[Dimension]
            return " * ".join(str(BaseUnits[d]) for d in dimension)

    @classmethod
    def get_contexts(cls, dimension: DimensionLike) -> list[Context]:
        if isinstance(dimension, str):
            return [BaseContexts[dimension]] if dimension in BaseContexts else []
        contexts = {BaseContexts[d] for d in dimension if d in BaseContexts}
        return list(contexts)

    @classmethod
    def from_pint(
        cls: "type[Self]",
        value: PlainQuantity[ValueType],
        dimension: DimensionLike,
    ) -> "Self":
        """Initialize a Quantity from a ureg.Quantity."""
        if isinstance(value, ureg.Quantity):
            return cls(value, dimension)
        elif isinstance(value, PlainQuantity):
            raise TypeError(
                "Only use pint quantities genereated by pairinteraction.ureg and not by a different pint.UnitRegistry."
            )
        else:
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
        value: Union[PlainQuantity[ValueType], ValueType],
        unit: Optional[str],
        dimension: DimensionLike,
    ) -> "Self":
        if unit is None:
            if np.isscalar(value) and np.all(value == 0):
                return cls.from_base_unit(value, dimension)
            return cls.from_pint(value, dimension)
        return cls.from_unit(value, unit, dimension)

    def to_pint(self) -> PlainQuantity[ValueType]:
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
        return self._quantity.to(unit, *contexts).magnitude

    def to_base_unit(self) -> ValueType:
        """Return the value of the quantity in the base unit."""
        return self.to_pint().to_base_units().magnitude

    def to_pint_or_unit(self, unit: Optional[str]) -> Union[ValueType, PlainQuantity[ValueType]]:
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
