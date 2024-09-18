"""Pydantic model for combined system with interactions."""

from typing import Optional, Union

from pydantic import BaseModel, ConfigDict, Field

from pairinteraction.model.states.atom import ModelStateAtomMQDT, ModelStateAtomSQDT
from pairinteraction.model.states.classical_light import ModelStateClassicalLight
from pairinteraction.model.types.parameter import (
    ParameterConstant,
    ParameterList,
    ParameterRange,
)
from pairinteraction.model.types.simple_types import ConstituentString, HalfInt, PositiveZero, Symmetry

UnionModelStates = Union[
    ModelStateAtomSQDT[int],
    ModelStateAtomSQDT[HalfInt],
    ModelStateAtomMQDT[int],
    ModelStateAtomMQDT[HalfInt],
    ModelStateClassicalLight,
]
UnionParameterInt = Union[ParameterConstant[int], ParameterList[int], ParameterRange[int]]
UnionParameterFloat = Union[ParameterConstant[float], ParameterList[float], ParameterRange[float]]
UnionParameterSymmetry = Union[ParameterConstant[Symmetry], ParameterList[Symmetry]]


class ModelInteractions(BaseModel):
    """Model corresponding to SystemWithInteractions."""

    model_config = ConfigDict(extra="forbid", frozen=False, validate_assignment=True)

    use_delta_energy_after_fields: bool = True

    min_energy: Optional[float] = None
    max_energy: Optional[float] = None
    delta_energy: Optional[PositiveZero[float]] = None

    min_energy_after_diagonalization: Optional[float] = None
    max_energy_after_diagonalization: Optional[float] = None
    delta_energy_after_diagonalization: Optional[PositiveZero[float]] = None

    conserved_total_m: Optional[UnionParameterInt] = None
    conserved_parity_under_inversion: UnionParameterSymmetry = Field(None, validate_default=True)
    conserved_parity_under_reflection: UnionParameterSymmetry = Field(None, validate_default=True)
    conserved_parity_under_permutation: UnionParameterSymmetry = Field(None, validate_default=True)

    distance: Optional[UnionParameterFloat] = None
    angle: UnionParameterFloat = Field(0, validate_default=True)

    combined_states_of_interest: list[dict[ConstituentString, Union[int, UnionModelStates]]] = []
