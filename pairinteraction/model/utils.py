"""Some extra utils, like special fields and validators for pydantic models,
that are going to be reused in different models."""

from pydantic import Field


def ExtraField(**kwargs):
    """Field with default=None and exclude=True."""
    return Field(None, exclude=True, **kwargs)


def LeftToRightField(*args, **kwargs):
    """Field with union_mode="left_to_right"."""
    return Field(*args, **kwargs, union_mode="left_to_right")
