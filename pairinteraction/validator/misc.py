"""Some more types, fields and other various useful stuff for the pydantic models."""

from typing import Union

from pydantic import Field, confloat, conint, conlist, constr
from typing_extensions import Literal

# TODO new python: replace conint, confloat, constr, conlist by Annotated[type, Field(...)]

PositiveFloat = confloat(ge=0)
PositiveInt = conint(ge=0)
HalfInt = confloat(multiple_of=0.5)
PositiveHalfInt = confloat(multiple_of=0.5, ge=0)

ConstituentString = Literal["atom1", "atom2", "classical_light1", "classical_light2"]
SpeciesString = constr(pattern=r"^[A-Z][a-z]?[0-9]{0,3}(_singlet|_triplet|_mqdt)?$")

Vector = conlist(float, min_length=3, max_length=3)

Symmetry = Literal[None, 1, -1]


def ExtraField(**kwargs):
    """Field with default=None and exclude=True."""
    return Field(None, exclude=True, **kwargs)


def LeftToRightField(*args, **kwargs):
    """Field with union_mode="left_to_right"."""
    return Field(*args, **kwargs, union_mode="left_to_right")


PossibleParameterTypes = Union[None, int, float]  # TODO add complex? for now this would throw a pydantic error
# TODO new python: needed for older python versions, for newer you can directly use PossibleParameterTypes
PossibleParameterTypesAsTuple = (type(None), int, float)

# TODO NEW: remove this
NEW_TO_OLD_SPECIES = {
    "Sr88_singlet": "Sr1",
    "Sr88_triplet": "Sr3",
}
OLD_TO_NEW_SPECIES = {v: k for k, v in NEW_TO_OLD_SPECIES.items()}
