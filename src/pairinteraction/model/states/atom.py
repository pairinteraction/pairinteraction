"""Pydantic models for states of atoms."""

from abc import ABC
from typing import Generic, Optional, TypeVar

from pydantic import Field, ValidationInfo, field_validator

from pairinteraction.model.states.base import BaseModelState
from pairinteraction.model.types.simple_types import (
    HalfInt,
    Positive,
    PositiveZero,
    QnTypes,
    SpeciesString,
    SpeciesStringMQDT,
    SpeciesStringSQDT,
)

SpinType = TypeVar("SpinType", int, HalfInt)
FType = TypeVar("FType", int, HalfInt)


class BaseModelStateAtom(BaseModelState, ABC):
    """Base model representing a state of a generic atom."""

    species: SpeciesString

    energy: QnTypes[float] = None


class ModelStateAtomSQDT(BaseModelStateAtom, Generic[SpinType]):
    """Model representing a state of a SQDT atom with either a integer or half integer spin
    (i.e. SpinType is int or HalfInt)."""

    species: SpeciesStringSQDT

    n: QnTypes[Positive[int]] = None
    nu: QnTypes[Positive[float]] = None
    l: QnTypes[PositiveZero[int]] = None
    s: QnTypes[PositiveZero[SpinType]] = Field(None, validate_default=True)
    j: QnTypes[PositiveZero[SpinType]] = None
    m: QnTypes[SpinType] = None

    @field_validator("s", mode="before")
    @classmethod
    def set_spin(cls, s: Optional[SpinType], info: ValidationInfo) -> SpinType:
        """If the spin is not provided, set it according to the species ending."""
        species = info.data.get("species")
        if species.endswith("singlet"):
            assert s is None or s == 0, "Spin must be 0 or species must end with '_singlet'."
            return 0
        elif species.endswith("triplet"):
            assert s is None or s == 1, "Spin must be 1 or species must end with '_triplet'."
            return 1
        else:
            assert s is None or s == 0.5, "Spin must be 0.5 or species must not end with '_singlet'/'_triplet'."
            return 0.5

    @property
    def f(self) -> QnTypes[PositiveZero[SpinType]]:
        """Return the total momentum quantum number, for SQDT states this is simply j."""
        return self.j


class ModelStateAtomMQDT(BaseModelStateAtom, Generic[FType]):
    """Model representing a state of a MQDT atom with either a integer or half integer total momentum f
    (i.e. FType is int or HalfInt).
    """

    species: SpeciesStringMQDT

    n: QnTypes[Positive[float]] = None
    nu: QnTypes[Positive[float]] = None
    l: QnTypes[PositiveZero[float]] = None
    s: QnTypes[PositiveZero[float]] = None
    j: QnTypes[PositiveZero[float]] = None
    f: QnTypes[PositiveZero[FType]] = None
    m: QnTypes[FType] = None
