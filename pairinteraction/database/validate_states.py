from typing import List

from pairinteraction import pireal
from pairinteraction.model.states import BaseModelState, ModelStateAtomSQDT


def validate_states_of_interest(states_of_interest: List[BaseModelState]) -> List[BaseModelState]:
    """Validate the states of interest and get the undefined quantum numbers and the energy from the database."""
    new_states_of_interest = []
    for state in states_of_interest:
        if not isinstance(state, ModelStateAtomSQDT):
            raise NotImplementedError("TODO.")
        new_states_of_interest += validate_sqdt_state(state)

    return new_states_of_interest


def validate_sqdt_state(state: ModelStateAtomSQDT) -> List[ModelStateAtomSQDT]:
    """Validate the SQDT state.

    TODO remove this with proper database implementation.
    """
    new_states = [state]

    for qn_name in ["n", "l", "s", "j", "m"]:
        qn = getattr(state, qn_name)
        if qn is None:
            raise ValueError(f"Quantum number {qn_name} must be provided.")
        elif isinstance(qn, list):
            old_states = new_states
            new_states = []
            for s in old_states:
                for qn_val in qn:
                    new_state = s.copy()
                    setattr(new_state, qn_name, qn_val)
                    new_states.append(new_state)

    for state in new_states:
        state.energy = pireal.StateOne(state.species, state.n, state.l, state.j, state.m).getEnergy()

    return new_states
