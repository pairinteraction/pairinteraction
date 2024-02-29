"""Pydantic models for parameters."""


from typing import Dict, List, Optional, Union

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, model_serializer, model_validator

from pairinteraction.validator.misc import ExtraField, PossibleParameterTypes


class BaseParameter(BaseModel):
    """Pydantic base model for a parameter.

    All fields, that should be able to be looped over should be of this type.
    """

    # TODO: is there a option to set frozen to True after all model_validators have been run?
    model_config = ConfigDict(extra="forbid", frozen=False)

    def get_value(self, step: Optional[int] = None) -> PossibleParameterTypes:
        raise NotImplementedError("This method should be implemented in the subclasses!")


class ParameterSimple(BaseParameter):
    """Pydantic model for a simple value as parameter."""

    value: PossibleParameterTypes = None

    @property
    def list(self) -> List[PossibleParameterTypes]:
        """Return the value as a list."""
        return [self.value]

    def get_value(self, step: Optional[int] = None) -> PossibleParameterTypes:
        """If the parameter is a ParameterSimple simply return the value, no matter the step."""
        return self.value

    @model_serializer(mode="wrap")
    def serialize_as_value(self, handler) -> PossibleParameterTypes:
        """Serialize the parameter as a simple value."""
        # TODO not sure why this is needed, see https://github.com/pydantic/pydantic/discussions/8541
        if isinstance(self, ParameterSimple):
            return self.value
        return handler(self)


class ParameterRange(BaseParameter):
    """Pydantic model for a parameter range."""

    start: Optional[float] = ExtraField()
    stop: Optional[float] = ExtraField()
    steps: Optional[int] = ExtraField()

    list: Optional[List[PossibleParameterTypes]] = None

    @model_validator(mode="after")
    def use_start_stop_steps(self) -> "ParameterRange":
        """If start, stop and steps is given, create the corresponding list using numpy.linspace."""
        start_stop_steps = [self.start, self.stop, self.steps]

        if self.list is not None:
            if self.start is not None or self.stop is not None:
                raise ValueError("If list is provided, start, stop and steps cannot be provided!")
            if self.steps is not None and self.steps != len(self.list):
                raise ValueError("If list is provided, steps cannot be provided (or must be exactle len(list)!")
            self.steps = len(self.list)
            return self

        if any(x is None for x in start_stop_steps):
            raise ValueError("Either list or (start, stop and steps) must be provided!")

        self.list = list(np.linspace(*start_stop_steps))
        self.start, self.stop = None, None  # TODO do i want to store those for nicer serialization?
        return self

    def get_value(self, step: Optional[int] = None) -> PossibleParameterTypes:
        """Get the value at a specific step."""
        if step is None:
            raise ValueError("ParameterRange: to get a step value a step must be provided!")
        return self.list[step]


UnionParameter = Union[ParameterSimple, ParameterRange]


class ParameterRangeOptions(BaseModel):
    """Pydantic model for collecting all parameter ranges."""

    model_config = ConfigDict(extra="forbid", frozen=False)

    steps: int = Field(default=None, ge=1)
    parameters: Dict[str, ParameterRange] = {}

    @model_validator(mode="after")
    def check_steps(self) -> "ParameterRangeOptions":
        """Check if steps for all parameters are the same and set self.steps if not provided."""
        if len(self.parameters) == 0:
            if not (self.steps is None or self.steps == 1):
                raise ValueError("If no parameters are given, steps must be None or 1!")
            self.steps = 1
            return self

        if self.steps is None:
            self.steps = next(iter(self.parameters.values())).steps

        for k, v in self.parameters.items():
            if self.steps != v.steps:
                raise ValueError(f"steps for parameter {k} is {v.steps}, but should be {self.steps}!")
        return self
