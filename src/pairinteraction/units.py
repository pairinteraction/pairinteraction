from collections.abc import Collection
from typing import TYPE_CHECKING, Any, Generic, Literal, TypeVar, Union

import numpy as np
import scipy.sparse
from pint import UnitRegistry
from pint.facets.plain import PlainQuantity

if TYPE_CHECKING:
    from pint.facets.plain import PlainUnit
    from scipy.sparse import csr_matrix

    Array = np.ndarray[Any, Any]
    SelfQuantity_t = TypeVar("SelfQuantity_t", bound="QuantityAbstract[Any]")

ureg = UnitRegistry(system="atomic")

Dimension = Literal["electric_field", "magnetic_field", "distance", "energy"]
BaseUnits: dict[Dimension, "PlainUnit"] = {
    "electric_field": ureg.Quantity(1, "V/cm").to_base_units().units,
    "magnetic_field": ureg.Quantity(1, "G").to_base_units().units,
    "distance": ureg.Quantity(1, "micrometer").to_base_units().units,
    "energy": ureg.Unit("hartree"),
}
ValueType = TypeVar("ValueType", bound=Union[float, "Array", "csr_matrix"])


class QuantityAbstract(Generic[ValueType]):
    def __init__(self, value: Union[ValueType, PlainQuantity[ValueType]], unit: Union[str, "PlainUnit"] = "pint"):
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
    def from_base(
        cls: "type[SelfQuantity_t]", value: Union[ValueType, PlainQuantity[ValueType]], dimension: Dimension
    ) -> "SelfQuantity_t":
        return cls(value, BaseUnits[dimension])

    def to_base(self, dimension: Dimension) -> ValueType:
        return self.to_unit(BaseUnits[dimension])  # type: ignore

    def to_unit(self, unit: Union[str, "PlainUnit"] = "pint") -> Union[ValueType, PlainQuantity[ValueType]]:
        if unit == "pint":
            return self.quantity
        return self.quantity.to(unit, "spectroscopy").magnitude  # type: ignore


class QuantityScalar(QuantityAbstract[float]):
    def __init__(self, value: Union[float, PlainQuantity[float]], unit: Union[str, "PlainUnit"] = "pint"):
        self._is_zero_without_unit = False
        if unit == "pint" and np.isscalar(value) and value == 0:
            value = ureg.Quantity(0)
            self._is_zero_without_unit = True
        super().__init__(value, unit)
        magnitude = self.quantity.magnitude
        if not np.isscalar(magnitude):
            raise TypeError(f"value must be a scalar, not {type(magnitude)}")

    def to_unit(self, unit: Union[str, "PlainUnit"] = "pint") -> Union[float, PlainQuantity[float]]:
        if self._is_zero_without_unit:
            if unit == "pint":
                raise ValueError("value 0 without a unit given cannot be converted to a pint.Quantity")
            return 0
        return super().to_unit(unit)


class QuantityArray(QuantityAbstract["Array"]):
    def __init__(self, value: Union["Array", PlainQuantity["Array"]], unit: Union[str, "PlainUnit"] = "pint"):
        super().__init__(value, unit)
        magnitude = self.quantity.magnitude
        if not isinstance(magnitude, Collection):
            raise TypeError(f"value must be an Array (or a Collection), not {type(magnitude)}")
        if not all(np.isscalar(v) for v in magnitude):
            raise TypeError(
                f"values must be a a list of scalars, not {type(magnitude)=}, {[type(v) for v in magnitude]=}"
            )


class QuantitySparse(QuantityAbstract["csr_matrix"]):
    def __init__(self, value: Union["csr_matrix", PlainQuantity["csr_matrix"]], unit: Union[str, "PlainUnit"] = "pint"):
        super().__init__(value, unit)
        magnitude = self.quantity.magnitude
        if not isinstance(magnitude, scipy.sparse.csr_matrix):
            raise TypeError(f"value must be a scipy.sparse.csr_matrix, not {type(magnitude)}")
