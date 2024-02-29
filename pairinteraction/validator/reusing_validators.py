"""Various validators, that are used multiple times in different models."""


from typing import Union

from pairinteraction.validator.misc import PossibleParameterTypesAsTuple


def one_use_delta_and_soi(self, attr: str, use_combined=False) -> None:
    """Check if delta_attr and (combined) states_of_interest are given and convert to min/max_attr for one attribute."""
    # TODO make this a field_validator instead of a model_validator
    delta = getattr(self, f"delta_{attr}")
    if delta is None:
        return

    min_max_given = [getattr(self, k) is not None for k in [f"min_{attr}", f"max_{attr}"]]
    if any(min_max_given):
        raise ValueError(f"delta_{attr} and (min_{attr} or max_{attr}) are given, use them exclusively.")

    states_of_interest = self.combined_states_of_interest if use_combined else self.states_of_interest
    if len(states_of_interest) == 0:
        raise ValueError(
            f"delta_{attr} given, but no states of interest are given. "
            "Also provide (combined_)states_of_interest or use min_{attr} and max_{attr} instead."
        )

    if use_combined and self.use_delta_energy_after_fields in [None, True]:
        self.use_delta_energy_after_fields = True
        return

    state_attr = "energy" if attr == "energy_after_diagonalization" else attr
    values = [getattr(state, state_attr) for state in states_of_interest]
    min_attr, max_attr = min(values) - delta, max(values) + delta

    if attr == "n":
        min_attr = max(min_attr, 1)
    elif attr == "l":
        min_attr = max(min_attr, 0)
    elif attr == "j" and min_attr < 0:
        min_attr = min_attr % 1

    setattr(self, f"min_{attr}", min_attr)
    setattr(self, f"max_{attr}", max_attr)


def use_parameter_if_float(value: Union[float, dict]) -> dict:
    """If value is already the parameter return a a simple dictionary {"value": value},
    so it can be parsed as a ParameterSimple.
    """
    if isinstance(value, PossibleParameterTypesAsTuple):
        return {"value": value}
    elif isinstance(value, list):
        return {"list": value}
    return value
