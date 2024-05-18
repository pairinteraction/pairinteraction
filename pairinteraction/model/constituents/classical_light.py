"""Pydantic model for the classical light."""

from typing import List, Optional

from pydantic import (
    ValidationInfo,
    field_validator,
)

from pairinteraction.model.constituents.base import BaseModelConstituent
from pairinteraction.model.states import ModelStateClassicalLight


class ModelClassicalLight(BaseModelConstituent):
    """Model representing classical light."""

    # TODO update this class

    frequency: float
    pointing_vector: List[float]
    intensity_per_polarization: List[float]

    min_n_pi: Optional[int] = None
    max_n_pi: Optional[int] = None
    min_n_plus: Optional[int] = None
    max_n_plus: Optional[int] = None
    min_n_minus: Optional[int] = None
    max_n_minus: Optional[int] = None
    min_n_total: Optional[int] = None
    max_n_total: Optional[int] = None

    min_energy_after_diagonalization: Optional[float] = None
    max_energy_after_diagonalization: Optional[float] = None

    states_of_interest: List[ModelStateClassicalLight] = []

    @field_validator("states_of_interest", mode="before")
    @classmethod
    def set_soi_attributes(cls, states_of_interest: List, info: ValidationInfo) -> List:
        """If states_of_interest are provided, automatically set (and check) the attributes
        (frequency, pointing_vector and intensity_per_polarization) of all the states.

        This allows to define the attributes only once and not for each state separately.
        """
        for soi in states_of_interest:
            for attr in ["frequency", "pointing_vector", "intensity_per_polarization"]:
                value = info.data.get(attr)
                assert (
                    soi.setdefault(attr, value) == value
                ), f"{attr} of states_of_interest must be the same as for the constituent"
        return states_of_interest
