from collections.abc import Collection
from typing import TYPE_CHECKING, Literal, Union

import numpy as np
from pint import UnitRegistry
from pint.facets.plain import PlainQuantity

if TYPE_CHECKING:
    from numpy.typing import ArrayLike
    from pint.facets.plain import PlainUnit

ureg = UnitRegistry(system="atomic")

Dimension = Literal["electric_field", "magnetic_field", "distance", "energy"]
BaseUnits: dict[Dimension, "PlainUnit"] = {
    "electric_field": ureg.Quantity(1, "V/cm").to_base_units().units,
    "magnetic_field": ureg.Quantity(1, "G").to_base_units().units,
    "distance": ureg.Quantity(1, "micrometer").to_base_units().units,
    "energy": ureg.Unit("hartree"),
}


class Qty:
    def __init__(self, value: Union[float, "PlainQuantity"], unit: Union[str, "PlainUnit"] = "pint"):
        if isinstance(value, ureg.Quantity):
            if unit != "pint":
                raise ValueError("unit must be 'pint' if value is a pairinteraction.ureg.Quantity")
            self.quantity = value
        elif isinstance(value, PlainQuantity):
            raise TypeError(
                "Only use pint quantities genereated by pairinteraction.ureg and not by a different pint.UnitRegistry."
            )
        elif np.isscalar(value):
            if unit == "pint":
                raise ValueError(
                    "unit must be a string specifiying the unit of the value if value is a scalar (e.g. unit='GHz')"
                )
            self.quantity = ureg.Quantity(value, unit)
        else:
            raise TypeError("value must be a scalar or a pairinteraction.ureg.Quantity")

    @classmethod
    def from_base(cls, value: float, dimension: Dimension) -> "Qty":
        return cls(value, BaseUnits[dimension])

    def to_base(self, dimension: Dimension) -> float:
        return self.quantity.to(BaseUnits[dimension], "spectroscopy").magnitude

    def to_unit(self, unit: Union[str, "PlainUnit"] = "pint") -> Union[float, "PlainQuantity[float]"]:
        if unit == "pint":
            return self.quantity
        return self.quantity.to(unit, "spectroscopy").magnitude


class Qties:
    def __init__(
        self, values: Union["PlainQuantity", Collection[float], "ArrayLike"], unit: Union[str, "PlainUnit"] = "pint"
    ):
        if isinstance(values, ureg.Quantity):
            if unit != "pint":
                raise ValueError("unit must be 'pint' if values is a pairinteraction.ureg.Quantity")
            self.quantities = values
        elif isinstance(values, PlainQuantity):
            raise TypeError(
                "Only use pint quantities genereated by pairinteraction.ureg and not by a different pint.UnitRegistry."
            )
        elif isinstance(values, Collection):
            if not all(np.isscalar(v) for v in values):
                raise TypeError("values must be a a list of scalars")
            if unit == "pint":
                raise ValueError(
                    "if values is a list of scalars, unit must be a string specifiying the unit (e.g. unit='GHz')"
                )
            self.quantities = ureg.Quantity(values, unit)
        else:
            raise TypeError(
                f"values must be a list of scalars or a pairinteraction.ureg.Quantity, but were of type {type(values)}"
            )

    @classmethod
    def from_base(cls, values: Union[Collection[float], "ArrayLike"], dimension: Dimension) -> "Qties":
        return cls(values, BaseUnits[dimension])

    def to_base(self, dimension: Dimension) -> "ArrayLike":
        return self.quantities.to(BaseUnits[dimension], "spectroscopy").magnitude

    def to_unit(self, unit: Union[str, "PlainUnit"]) -> Union["ArrayLike", "PlainQuantity"]:
        if unit == "pint":
            return self.quantities
        return self.quantities.to(unit, "spectroscopy").magnitude
