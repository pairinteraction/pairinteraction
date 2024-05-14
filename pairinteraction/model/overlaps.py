"""Pydantic model for defining what kind of overlaps to calculate."""


from typing import Dict, List, Optional, Union

from pydantic import BaseModel, ConfigDict

from pairinteraction.model.states import UnionModelStates
from pairinteraction.model.types import ConstituentString


class ModelOverlaps(BaseModel):
    """Pydantic model for the overlaps between states."""

    model_config = ConfigDict(extra="forbid", frozen=False)

    combined_states_of_interest: List[Dict[ConstituentString, Union[int, UnionModelStates]]] = []
    transformation: Dict[ConstituentString, Optional[str]] = {}
