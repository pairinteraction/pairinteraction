from typing import List

from pairinteraction import pireal
from pairinteraction.model.constituents.atom import BaseModelAtom
from pairinteraction.model.constituents.base import BaseModelConstituent
from pairinteraction.model.simulation import ModelSimulation
from pairinteraction.model.states.atom import ModelStateAtomSQDT


def preprocess_model_simulation(model: ModelSimulation) -> None:
    """Preprocess the model simulation."""
    for constituent in model.unique_constituents.values():
        preprocess_states(constituent)
        preprocess_constituent_restrictions(constituent)

    if model.interactions is not None:
        preprocess_interaction_restrictions(model)


def preprocess_constituent_restrictions(constituent: BaseModelConstituent) -> None:
    """If states_of_interest is provided together with delta_attr instead of min/max_attr
    convert the delta_attr to min/max_attr.
    """
    for attr in ["n", "nu", "l", "s", "j", "f", "m", "energy"]:
        preprocess_one_atom_restrictions(constituent, attr)


def preprocess_one_atom_restrictions(atom: BaseModelAtom, attr: str) -> None:
    """Convert delta_attr to min/max_attr for a single atom and one attribute."""
    delta = getattr(atom, f"delta_{attr}", None)
    if delta is None:
        return

    min_, max_ = getattr(atom, f"min_{attr}"), getattr(atom, f"max_{attr}")
    if (min_ is None or max_ is None) and min_ != max_:
        raise ValueError(
            f"For atom delta_{attr} given and (min_{attr} xor max_{attr}) is given. " "This behaviour is not defined."
        )

    states_of_interest = atom.states_of_interest
    if len(states_of_interest) == 0:
        raise ValueError(
            f"delta_{attr} given, but no states of interest are given. "
            f"Also provide states_of_interest or use min_{attr} and max_{attr} instead."
        )

    values = [getattr(state, attr, None) for state in states_of_interest]
    if any(v is None for v in values):
        raise ValueError(
            f"delta_{attr} given, but not all states of interest have a {attr}."
            "Consider calling preprocess_states before preprocess_constituent_restrictions."
        )
    min_attr, max_attr = min(values) - delta, max(values) + delta

    if attr == "n":
        min_attr = max(min_attr, 1)
    elif attr == "l":
        min_attr = max(min_attr, 0)
    elif attr == "j":
        state0 = states_of_interest[0]
        if isinstance(state0, ModelStateAtomSQDT):
            min_attr = max(min_attr, state0.s % 1)
        else:
            min_attr = max(min_attr, 0)

    for minmax, new in zip(["min", "max"], [min_attr, max_attr]):
        old = getattr(atom, f"{minmax}_{attr}")
        if old is None:
            setattr(atom, f"{minmax}_{attr}", new)
        elif old != new:
            raise ValueError(f"delta_{attr} given and {minmax}_{attr} is already set but they are not compatible.")


def preprocess_interaction_restrictions(model: ModelSimulation) -> ModelSimulation:
    # TODO how to handle this (also in the context of applying delta after diagonalizing single atoms)
    pass


def preprocess_states(constituent: BaseModelConstituent) -> None:
    """Validate the states of interest and get the undefined quantum numbers and the energy from the database."""
    new_states_of_interest = []
    for state in constituent.states_of_interest:
        if not isinstance(state, ModelStateAtomSQDT):
            raise NotImplementedError("TODO.")
        new_states_of_interest += preprocess_sqdt_state(state)
    constituent.states_of_interest = new_states_of_interest


def preprocess_sqdt_state(state: ModelStateAtomSQDT) -> List[ModelStateAtomSQDT]:
    """Validate the SQDT state.

    TODO NEW remove this with proper database implementation.
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
