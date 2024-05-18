from abc import ABC

from pydantic import BaseModel, ConfigDict


class BaseModelState(BaseModel, ABC):
    """Base model representing a state of a generic constituent."""

    model_config = ConfigDict(extra="forbid", frozen=False, validate_assignment=True)
