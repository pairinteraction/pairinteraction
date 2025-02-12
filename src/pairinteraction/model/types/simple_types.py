"""Define some useful types and unions."""

from typing import Annotated, Any, Literal, TypeVar, Union

from pydantic import Field, GetCoreSchemaHandler
from pydantic_core import core_schema

T = TypeVar("T")


ConstituentString = Literal["atom1", "atom2"]

SpeciesString = Annotated[str, Field(pattern=r"^[A-Z][a-z]?[0-9]{0,3}(_singlet|_triplet|_mqdt)?$")]
SpeciesStringSQDT = Annotated[str, Field(pattern=r"^[A-Z][a-z]?[0-9]{0,3}(_singlet|_triplet)?$")]
SpeciesStringMQDT = Annotated[str, Field(pattern=r"^[A-Z][a-z]?[0-9]{0,3}_mqdt$")]

Vector = Annotated[list[float], Field(min_length=3, max_length=3)]

Symmetry = Literal[None, 1, -1]


class HalfInt(float):
    """Class for defining half integer."""

    def __init__(self, value: float) -> None:
        """Check that value is a half integer (i.e. value % 1 == 0.5)."""
        if value % 1 != 0.5:
            raise ValueError(f"{value} is not a half integer")

    @classmethod
    def __get_pydantic_core_schema__(cls, source: type[Any], handler: GetCoreSchemaHandler) -> core_schema.CoreSchema:
        assert source is HalfInt
        return core_schema.no_info_after_validator_function(
            cls,
            core_schema.float_schema(),
        )


Positive = Annotated[T, Field(gt=0)]
PositiveZero = Annotated[T, Field(ge=0)]

QnTypes = Union[T, list[T], Literal["any"], None]  # TODO allow range from dict
