"""Some extra utils, like special fields and validators for pydantic models,
that are going to be reused in different models."""


import collections.abc
from typing import Union, get_args

from pydantic import Field

from pairinteraction.model.types import PossibleParameterTypes


def ExtraField(**kwargs):
    """Field with default=None and exclude=True."""
    return Field(None, exclude=True, **kwargs)


def LeftToRightField(*args, **kwargs):
    """Field with union_mode="left_to_right"."""
    return Field(*args, **kwargs, union_mode="left_to_right")


def validate_parameter(p: Union[PossibleParameterTypes, dict]) -> dict:
    """If p is a simple value return {"value": value} or if it is a list return {"list": list},
    so it can be parsed by the pydantic Parameter classes derived from BaseParameter.
    """
    if isinstance(p, get_args(PossibleParameterTypes)):
        return {"value": p}
    elif isinstance(p, collections.abc.Sequence) and not isinstance(p, str):
        return {"list": p}
    return p
