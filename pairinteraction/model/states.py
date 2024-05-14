"""Pydantic models for states of constituents."""

from functools import cached_property
from typing import List, Optional, Union

from pydantic import BaseModel, ConfigDict, field_validator

from pairinteraction import pireal
from pairinteraction.model.types import (
    HalfInt,
    PositiveFloat,
    PositiveHalfInt,
    PositiveInt,
    SpeciesString,
)


def AllowedQn(type_):
    """Return a Union of all allowed types for a quantum number.

    We allow the quantum number to be a single value, a list of values or None.
    """
    return Union[type_, List[type_], None]


class BaseModelState(BaseModel):
    """Pydantic base model representing a state of a generic constituent."""

    model_config = ConfigDict(extra="forbid", frozen=True)


class BaseModelStateAtom(BaseModelState):
    """Pydantic base model representing a state of a generic atom (simple or mqdt)."""

    # TODO NEW: link state to database and ID (where to get the database_directory from?)
    # then energy, s, ... can simple be looked up from the database
    species: SpeciesString

    @cached_property
    def energy(self) -> float:
        """Get the energy of the state from the database."""
        # TODO NEW: get energy from database (also remove the methods in the derived classes)
        raise NotImplementedError("Get energy from database not implemented yet.")


class ModelStateAtomSimple(BaseModelStateAtom):
    """Pydantic model representing a state of a simple (=no mqdt) atom."""

    n: AllowedQn(PositiveInt) = None
    l: AllowedQn(PositiveInt) = None
    j: AllowedQn(Union[PositiveInt, PositiveHalfInt]) = None
    m: AllowedQn(Union[int, HalfInt]) = None

    @cached_property
    def s(self) -> PositiveHalfInt:
        """Get the spin of the atom from the provided species."""
        species = self.species
        no_ending = not any(species.endswith(ending) for ending in ["singlet", "triplet", "mqdt"])

        if no_ending:
            return 1 / 2
        if species.endswith("singlet"):
            return 0
        if species.endswith("triplet"):
            return 1
        raise ValueError(f"Never reached: Species {species} is not a valid species for a simple atom.")

    # no computed_field -> excluded from model_dump
    @property
    def f(self) -> Optional[HalfInt]:
        """Return the total momentum quantum number, for simple (non MQDT) states this is simply self.j ."""
        return self.j

    @cached_property
    def energy(self) -> float:
        """Calculate the energy of the state."""
        m = self.m if self.m is not None else self.j
        if any(x is None or isinstance(x, list) for x in [self.n, self.l, self.j]):
            raise ValueError(
                "The quantum numbers n, l and j must be provided as single value"
                "(not as list or None) to calculate the energy."
            )
        return pireal.StateOne(self.species, self.n, self.l, self.j, m).getEnergy()

    @cached_property
    def state(self) -> pireal.StateOne:
        """Return the pairinteraction state object."""
        if any(x is None or isinstance(x, list) for x in [self.n, self.l, self.j]):
            raise ValueError(
                "The quantum numbers n, l, j and m must be provided as single value"
                "(not as list or None) to get the state."
            )
        return pireal.StateOne(self.species, self.n, self.l, self.j, self.m)

    @field_validator("species", mode="after")
    @classmethod
    def check_species(cls, species: SpeciesString) -> SpeciesString:
        """Check if the species is a valid simple atom."""
        if species.endswith("mqdt"):
            raise ValueError(f"Species {species} is not a valid species for a simple atom.")
        return species


class ModelStateAtomMQDT(BaseModelStateAtom):
    """Pydantic model representing a state of a atom using MQDT."""

    n: AllowedQn(PositiveFloat) = None
    nu: AllowedQn(PositiveFloat) = None
    l: AllowedQn(PositiveFloat) = None
    s: AllowedQn(PositiveFloat) = None
    j: AllowedQn(PositiveFloat) = None
    f: AllowedQn(PositiveHalfInt) = None
    m: AllowedQn(HalfInt) = None


class ModelStateClassicalLight(BaseModelState):
    """Pydantic model for the states of interest for ModelClassicalLight."""

    frequency: float
    pointing_vector: List[float]
    intensity_per_polarization: List[float]
    n_pi: Optional[int] = None
    n_plus: Optional[int] = None
    n_minus: Optional[int] = None
    n_total: Optional[int] = None

    @cached_property
    def energy(self) -> float:
        """Calculate the energy of the classical light state."""
        # TODO
        raise NotImplementedError("Energy of classical light state not implemented yet.")


UnionModelStates = Union[ModelStateAtomSimple, ModelStateAtomMQDT, ModelStateClassicalLight]
