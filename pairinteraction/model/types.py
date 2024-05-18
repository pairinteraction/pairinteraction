"""Define some useful types and unions."""

from typing import Any, List, Literal, Type, TypeVar, Union

from pydantic import Field, GetCoreSchemaHandler
from pydantic_core import core_schema
from typing_extensions import Annotated

T = TypeVar("T")


ConstituentString = Literal["atom1", "atom2", "classical_light1", "classical_light2"]

SpeciesString = Annotated[str, Field(pattern=r"^[A-Z][a-z]?[0-9]{0,3}(_singlet|_triplet|_mqdt)?$")]
SpeciesStringSQDT = Annotated[str, Field(pattern=r"^[A-Z][a-z]?[0-9]{0,3}(_singlet|_triplet)?$")]
SpeciesStringMQDT = Annotated[str, Field(pattern=r"^[A-Z][a-z]?[0-9]{0,3}_mqdt$")]

Vector = Annotated[List[float], Field(min_length=3, max_length=3)]

Symmetry = Literal[None, 1, -1]


class HalfIntAnnotation:
    """Class for defining pydantic schema for half integer."""

    def __get_pydantic_core_schema__(self, source: Type[Any], handler: GetCoreSchemaHandler) -> core_schema.CoreSchema:
        schema = handler(source)
        if schema["type"] != "float":
            raise TypeError("HalfIntAnnotation can only be applied to float")
        return core_schema.no_info_after_validator_function(
            self.validate,
            schema,
        )

    def validate(self, value: float) -> float:
        """Validate that value is a half integer (i.e. value % 1 == 0.5)."""
        if value % 1 != 0.5:
            raise ValueError(f"{value} is not a half integer")
        return value


HalfInt = Annotated[float, HalfIntAnnotation()]

Positive = Annotated[T, Field(gt=0)]
PositiveZero = Annotated[T, Field(ge=0)]

QnTypes = Union[T, List[T], Literal["any"], None]  # TODO allow range from dict


def ExcludedField(*args, **kwargs):
    """Field with exclude=True."""
    return Field(*args, exclude=True, **kwargs)


def ValidatedField(*args, **kwargs):
    """Field with validate_default=True."""
    return Field(*args, validate_default=True, **kwargs)
