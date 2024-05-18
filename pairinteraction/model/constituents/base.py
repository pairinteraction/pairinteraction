from abc import ABC

from pydantic import BaseModel, ConfigDict


class BaseModelConstituent(BaseModel, ABC):
    """Base model for the constituents.

    Each constituent corresponds to a System... .
    """

    model_config = ConfigDict(extra="forbid", frozen=False, validate_assignment=True)
