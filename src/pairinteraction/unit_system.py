from collections.abc import Sequence
from typing import TYPE_CHECKING, Literal, Union

import numpy as np
from numpy.typing import ArrayLike
from pint import UnitRegistry

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity, PlainUnit

ureg = UnitRegistry(system="atomic")


Dimension = Literal["electric_field", "magnetic_field", "distance", "energy"]
BaseUnits: dict[Dimension, "PlainUnit"] = {
    "electric_field": ureg.Quantity(1, "V/cm").to_base_units().units,
    "magnetic_field": ureg.Quantity(1, "G").to_base_units().units,
    "distance": ureg.Quantity(1, "micrometer").to_base_units().units,
    "energy": ureg.Unit("hartree"),
}


def convert_quantity_to_raw(
    quantity: Union["PlainQuantity", tuple[float, str]], unit_to: Union[str, "PlainUnit"]
) -> Union[float, np.ndarray]:
    """Convert a pint.Quantity object to a raw float value in the given unit."""
    if isinstance(quantity, (float, int, complex)) and quantity == 0:
        return 0
    if isinstance(quantity, Sequence) and len(quantity) == 2:
        quantity = ureg.Quantity(*quantity)
    if not isinstance(quantity, ureg.Quantity):
        raise TypeError("quantity must be a pairinteraction.ureg.Quantity or a tuple[float, str] of (value, unit)")
    return quantity.to(unit_to, "spectroscopy").magnitude


def convert_quantity_to_base(
    quantity: Union["PlainQuantity", tuple[float, str]], dimension: Dimension
) -> Union[float, np.ndarray]:
    """Convert a pint.Quantity object to a raw float value in the base unit of the given dimension."""
    if dimension not in BaseUnits:
        raise ValueError(f"dimension must be one of {BaseUnits.keys()}")
    return convert_quantity_to_raw(quantity, BaseUnits[dimension])


def convert_raw_to_base(
    value: Union[float, Sequence[float], ArrayLike], unit_from: Union[str, "PlainUnit"], dimension: Dimension
) -> Union[float, np.ndarray]:
    """Convert a raw float value with a given unit to the base unit of the given dimension."""
    if not isinstance(unit_from, (str, ureg.Unit)):
        raise TypeError("unit_from must be a string or a pairinteraction.ureg.Unit")
    quantity = ureg.Quantity(value, unit_from)
    return convert_quantity_to_base(quantity, dimension)


def convert_base_to_quantity(
    value: Union[float, Sequence[float], ArrayLike], dimension: Dimension, unit_to: Union[str, "PlainUnit", None]
) -> "PlainQuantity":
    """Convert a raw float value in the base unit to a pint.Quantity object in the given unit."""
    if dimension not in BaseUnits:
        raise ValueError(f"dimension must be one of {BaseUnits.keys()}")
    quantity = ureg.Quantity(value, BaseUnits[dimension])
    if unit_to is None:
        return quantity
    return quantity.to(unit_to, "spectroscopy")
