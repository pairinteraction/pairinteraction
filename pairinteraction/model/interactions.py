"""Pydantic model for combined system with interactions."""

from typing import Dict, List, Optional, Union

from pydantic import BaseModel, ConfigDict, Field

from pairinteraction.model.parameter import UnionParameterFloat, UnionParameterInt, UnionParameterSymmetry
from pairinteraction.model.states import UnionModelStates
from pairinteraction.model.types import ConstituentString


class ModelInteractions(BaseModel):
    """Pydantic model corresponding to SystemWithInteractions."""

    # TODO: is there a option to set frozen to True after all model_validators have been run?
    model_config = ConfigDict(extra="forbid", frozen=False, validate_assignment=True)

    min_energy: Optional[float] = None
    max_energy: Optional[float] = None
    min_energy_after_diagonalization: Optional[float] = None
    max_energy_after_diagonalization: Optional[float] = None

    conserved_total_m: Optional[UnionParameterInt] = None
    conserved_parity_under_inversion: UnionParameterSymmetry = Field(None, validate_default=True)
    conserved_parity_under_reflection: UnionParameterSymmetry = Field(None, validate_default=True)
    conserved_parity_under_permutation: UnionParameterSymmetry = Field(None, validate_default=True)

    distance: UnionParameterFloat = Field(None, validate_default=True)
    angle: UnionParameterFloat = Field(0, validate_default=True, description="Angle between the two states in degrees.")

    # Optional use delta_attr and combined_states_of_interest to define min_/max_attr
    # dont use BaseModelState, but UnionModelStates to enable automatic parsing (same for UnionModelConstituent)
    combined_states_of_interest: List[Dict[ConstituentString, Union[int, UnionModelStates]]] = []
    # TODO both delta_energy(...) should be ExtraField(), but how to handle use_delta_energy_after_fields?
    delta_energy: Optional[float] = None
    delta_energy_after_diagonalization: Optional[float] = None
    # TODO: change name of use_delta_energy_after_fields?
    use_delta_energy_after_fields: bool = Field(True)

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
