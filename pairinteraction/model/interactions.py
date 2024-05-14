"""Pydantic model for combined system with interactions."""

from typing import Dict, List, Optional, Union

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from pairinteraction.model.parameter import UnionParameter
from pairinteraction.model.states import UnionModelStates
from pairinteraction.model.types import ConstituentString
from pairinteraction.model.utils import (
    one_use_delta_and_soi,
    validate_parameter,
)


class ModelInteractions(BaseModel):
    """Pydantic model corresponding to SystemWithInteractions."""

    # TODO: is there a option to set frozen to True after all model_validators have been run?
    model_config = ConfigDict(extra="forbid", frozen=False, validate_assignment=True)

    min_energy: Optional[float] = None
    max_energy: Optional[float] = None
    min_energy_after_diagonalization: Optional[float] = None
    max_energy_after_diagonalization: Optional[float] = None

    conserved_total_m: Optional[int] = None
    conserved_parity_under_inversion: UnionParameter = Field(None, validate_default=True)
    conserved_parity_under_reflection: UnionParameter = Field(None, validate_default=True)
    conserved_parity_under_permutation: UnionParameter = Field(None, validate_default=True)

    distance: UnionParameter = Field(None, validate_default=True)
    angle: UnionParameter = Field(0, validate_default=True, description="Angle between the two states in degrees.")

    # Optional use delta_attr and combined_states_of_interest to define min_/max_attr
    # dont use BaseModelState, but UnionModelStates to enable automatic parsing (same for UnionModelConstituent)
    _used_delta: bool = False
    combined_states_of_interest: List[Dict[ConstituentString, Union[int, UnionModelStates]]] = []
    # TODO both delta_energy(...) should be ExtraField(), but how to handle use_delta_energy_after_fields?
    delta_energy: Optional[float] = None
    delta_energy_after_diagonalization: Optional[float] = None
    # TODO: change name of use_delta_energy_after_fields?
    use_delta_energy_after_fields: Optional[bool] = None

    validate_parameter = field_validator(
        "distance",
        "angle",
        "conserved_parity_under_inversion",
        "conserved_parity_under_reflection",
        "conserved_parity_under_permutation",
        mode="before",
    )(validate_parameter)

    @model_validator(mode="after")
    def use_delta_and_csoi(self) -> "ModelInteractions":
        """If combined_states_of_interest is provided together with delta_attr instead of min/max_attr
        convert delta_attr it to min/max_attr .
        """
        # FIXME: this is a ugly hack to avoid running use_delta again when setting
        # and thus revalidating the min/max values
        if self._used_delta:
            return self
        self._used_delta = True
        for attr in ["energy", "energy_after_diagonalization"]:
            one_use_delta_and_soi(self, attr, use_combined=True)
        return self

    # def validate_conserved_total_m(cls, v: float) -> float:
    #     # TODO allow for list of pair momenta?
    #     # check if all combined_states_of_interest have the same m
    #     # check angle/quantization_axis/fields if momentum can or cannot be conserved?
    #     raise NotImplementedError

    # @model_validator(mode="after")
    # def apply_angle(self) -> "ModelInteractions":
    #     if self.angle is None or self.angle == 0:
    #         self.angle = None
    #         return self

    #     if isinstance(self.angle, float):
    #         raise NotImplementedError
    #         # TODO once this is implemented, maybe use ExtraField() for angle?
    #         # (not sure because still the question how to handle range of angle?)

    #         self.angle = None
    #         return self

    #     raise NotImplementedError("How to handle list of angle?")
