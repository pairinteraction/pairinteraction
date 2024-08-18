"""Pydantic model for defining what kind of overlaps to calculate."""

from typing import Dict, List, Optional, Union

from pydantic import BaseModel, ConfigDict

from pairinteraction.model.states.atom import ModelStateAtomMQDT, ModelStateAtomSQDT
from pairinteraction.model.states.classical_light import ModelStateClassicalLight
from pairinteraction.model.types.simple_types import ConstituentString, HalfInt

UnionModelStates = Union[
    ModelStateAtomSQDT[int],
    ModelStateAtomSQDT[HalfInt],
    ModelStateAtomMQDT[int],
    ModelStateAtomMQDT[HalfInt],
    ModelStateClassicalLight,
]


class ModelOverlaps(BaseModel):
    """Model for calculating the overlaps between states."""

    model_config = ConfigDict(extra="forbid", frozen=False, validate_assignment=True)

    combined_states_of_interest: List[Dict[ConstituentString, Union[int, UnionModelStates]]] = []
    transformation: Dict[ConstituentString, Optional[str]] = {}
