"""Pydantic model for a state of classical light."""


from typing import List, Optional

from pairinteraction.model.states.base import BaseModelState


class ModelStateClassicalLight(BaseModelState):
    """Model representing a state of classical light."""

    # TODO update this class

    frequency: float
    pointing_vector: List[float]
    intensity_per_polarization: List[float]
    n_pi: Optional[int] = None
    n_plus: Optional[int] = None
    n_minus: Optional[int] = None
    n_total: Optional[int] = None
