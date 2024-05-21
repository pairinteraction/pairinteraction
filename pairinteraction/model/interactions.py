"""Pydantic model for combined system with interactions."""

from typing import Dict, List, Optional, Union

from pydantic import BaseModel, ConfigDict, Field

from pairinteraction.model.states import ModelStateAtomMQDT, ModelStateAtomSQDT, ModelStateClassicalLight
from pairinteraction.model.types import (
    ConstituentString,
    HalfInt,
    ParameterConstantFloat,
    ParameterConstantInt,
    ParameterConstantSymmetry,
    ParameterListFloat,
    ParameterListInt,
    ParameterListSymmetry,
    ParameterRangeFloat,
    ParameterRangeInt,
    PositiveZero,
)

UnionModelStates = Union[
    ModelStateAtomSQDT[int],
    ModelStateAtomSQDT[HalfInt],
    ModelStateAtomMQDT[int],
    ModelStateAtomMQDT[HalfInt],
    ModelStateClassicalLight,
]
UnionParameterInt = Union[ParameterConstantInt, ParameterListInt, ParameterRangeInt]
UnionParameterFloat = Union[ParameterConstantFloat, ParameterListFloat, ParameterRangeFloat]
UnionParameterSymmetry = Union[ParameterConstantSymmetry, ParameterListSymmetry]


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
    angle: UnionParameterFloat = ParameterConstantFloat(0)

    combined_states_of_interest: List[Dict[ConstituentString, Union[int, UnionModelStates]]] = []
