from collections.abc import Collection
from typing import TYPE_CHECKING, Any, Generic, Literal, Optional, TypeVar, Union

import numpy as np
import scipy.sparse
from pint import UnitRegistry
from pint.facets.plain import PlainQuantity

if TYPE_CHECKING:
    from pint.facets.plain import PlainUnit
    from scipy.sparse import csr_matrix
    from typing_extensions import Self

    Array = np.ndarray[Any, Any]

ureg = UnitRegistry(system="atomic")

Dimension = Literal[
    "ELECTRIC_FIELD",
    "MAGNETIC_FIELD",
    "DISTANCE",
    "ENERGY",
    "ELECTRIC_DIPOLE",
    "ELECTRIC_QUADRUPOLE",
    "ELECTRIC_QUADRUPOLE_ZERO",
    "ELECTRIC_OCTUPOLE",
    "MAGNETIC_DIPOLE",
    "ARBITRARY",
    "ZERO",
]
BaseUnits: dict[Dimension, "PlainUnit"] = {
    # ELECTRIC_FIELD: 1 V/cm = 1.94469038e-10 electron_mass * bohr / atomic_unit_of_time ** 3 / atomic_unit_of_current
    "ELECTRIC_FIELD": ureg.Quantity(1, "V/cm").to_base_units().units,
    # MAGNETIC_FIELD: 1 T = 4.25438216e-06 electron_mass / atomic_unit_of_time ** 2 / atomic_unit_of_current
    "MAGNETIC_FIELD": ureg.Quantity(1, "T").to_base_units().units,
    # DISTANCE: 1 mum = 18897.2612 bohr
    "DISTANCE": ureg.Quantity(1, "micrometer").to_base_units().units,
    # ENERGY: 1 hartree = 1 electron_mass * bohr ** 2 / atomic_unit_of_time ** 2
    "ENERGY": ureg.Unit("hartree"),
    # ELECTRIC_DIPOLE: 1 e * a0 = 1 atomic_unit_of_current * atomic_unit_of_time * bohr
    "ELECTRIC_DIPOLE": ureg.Quantity(1, "e * a0").to_base_units().units,
    # ELECTRIC_QUADRUPOLE: 1 e * a0^2 = 1 atomic_unit_of_current * atomic_unit_of_time * bohr ** 2
    "ELECTRIC_QUADRUPOLE": ureg.Quantity(1, "e * a0^2").to_base_units().units,
    # ELECTRIC_QUADRUPOLE_ZERO: 1 e * a0^2 = 1 atomic_unit_of_current * atomic_unit_of_time * bohr ** 2
    "ELECTRIC_QUADRUPOLE_ZERO": ureg.Quantity(1, "e * a0^2").to_base_units().units,
    # ELECTRIC_OCTUPOLE: 1 e * a0^3 = 1 atomic_unit_of_current * atomic_unit_of_time * bohr ** 3
    "ELECTRIC_OCTUPOLE": ureg.Quantity(1, "e * a0^3").to_base_units().units,
    # MAGNETIC_DIPOLE: 1 hbar e / m_e = 1 bohr ** 2 * atomic_unit_of_current
    "MAGNETIC_DIPOLE": ureg.Quantity(1, "hbar e / m_e").to_base_units().units,
    "ARBITRARY": ureg.Unit(""),
    "ZERO": ureg.Unit(""),
}
Context = Literal["spectroscopy", "Gaussian"]
BaseContexts: dict[Dimension, Context] = {
    "MAGNETIC_FIELD": "Gaussian",
    "ENERGY": "spectroscopy",
}
ValueType = TypeVar("ValueType", bound=Union[float, "Array", "csr_matrix"])


class QuantityAbstract(Generic[ValueType]):
    def __init__(
        self, value: Union[ValueType, PlainQuantity[ValueType]], unit: Union[str, "PlainUnit"] = "pint"
    ) -> None:
        if isinstance(value, ureg.Quantity):
            if unit != "pint":
                raise ValueError("unit must be 'pint' if value is a pairinteraction.ureg.Quantity")
            self.quantity = value
        elif isinstance(value, PlainQuantity):
            raise TypeError(
                "Only use pint quantities genereated by pairinteraction.ureg and not by a different pint.UnitRegistry."
            )
        elif unit != "pint":
            self.quantity = ureg.Quantity(value, unit)
        else:  # value is not a pint.Quantity and unit is "pint"
            raise ValueError("unit must be a string specifiying the unit of the value if value is not a pint.Quantity")

    @classmethod
    def from_base(cls: "type[Self]", value: Union[ValueType, PlainQuantity[ValueType]], dimension: Dimension) -> "Self":
        return cls(value, BaseUnits[dimension])

    def to_base(self, dimension: Dimension) -> ValueType:
        context = BaseContexts.get(dimension, None)
        return self.to_unit(BaseUnits[dimension], context)  # type: ignore

    def to_unit(
        self, unit: Union[str, "PlainUnit"] = "pint", context: Optional[Context] = None
    ) -> Union[ValueType, PlainQuantity[ValueType]]:
        if unit == "pint":
            return self.quantity
        if context is None:
            dimension: list[Dimension] = [
                d for d, u in BaseUnits.items() if u.dimensionality == self.quantity.dimensionality
            ]
            if len(dimension) > 0:
                context = BaseContexts.get(dimension[0], None)
        if context is None:
            return self.quantity.to(unit).magnitude
        return self.quantity.to(unit, context).magnitude


class QuantityScalar(QuantityAbstract[float]):
    def __init__(self, value: Union[float, PlainQuantity[float]], unit: Union[str, "PlainUnit"] = "pint") -> None:
        self._is_zero_without_unit = False
        if unit == "pint" and np.isscalar(value) and value == 0:
            value = ureg.Quantity(0)
            self._is_zero_without_unit = True
        super().__init__(value, unit)
        magnitude = self.quantity.magnitude
        if not np.isscalar(magnitude):
            raise TypeError(f"value must be a scalar, not {type(magnitude)}")

    def to_unit(
        self, unit: Union[str, "PlainUnit"] = "pint", context: Optional[Context] = None
    ) -> Union[float, PlainQuantity[float]]:
        if self._is_zero_without_unit:
            if unit == "pint":
                raise ValueError("value 0 without a unit given cannot be converted to a pint.Quantity")
            return 0
        return super().to_unit(unit, context)


class QuantityArray(QuantityAbstract["Array"]):
    def __init__(self, value: Union["Array", PlainQuantity["Array"]], unit: Union[str, "PlainUnit"] = "pint") -> None:
        super().__init__(value, unit)
        magnitude = self.quantity.magnitude
        if not isinstance(magnitude, Collection):
            raise TypeError(f"value must be an Array (or a Collection), not {type(magnitude)}")
        if not all(np.isscalar(v) for v in magnitude):
            raise TypeError(
                f"values must be a a list of scalars, not {type(magnitude)=}, {[type(v) for v in magnitude]=}"
            )


class QuantitySparse(QuantityAbstract["csr_matrix"]):
    def __init__(
        self, value: Union["csr_matrix", PlainQuantity["csr_matrix"]], unit: Union[str, "PlainUnit"] = "pint"
    ) -> None:
        super().__init__(value, unit)
        magnitude = self.quantity.magnitude
        if not isinstance(magnitude, scipy.sparse.csr_matrix):
            raise TypeError(f"value must be a scipy.sparse.csr_matrix, not {type(magnitude)}")
