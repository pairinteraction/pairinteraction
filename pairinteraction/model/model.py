"""Pydantic model for the complete simulation."""

from functools import cached_property
from typing import Dict, Optional, Union

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    FieldSerializationInfo,
    SerializerFunctionWrapHandler,
    ValidationInfo,
    field_serializer,
    field_validator,
    model_validator,
)

from pairinteraction.model.constituents import (
    BaseModelConstituent,
    ModelAtom,
    ModelClassicalLight,
)
from pairinteraction.model.interactions import ModelInteractions
from pairinteraction.model.numerics import ModelNumerics
from pairinteraction.model.overlaps import ModelOverlaps
from pairinteraction.model.parameter import ParameterRange, ParameterRangeOptions
from pairinteraction.model.types import ConstituentString


class ModelSimulation(BaseModel):
    """Pydantic model for the input of a pairinteraction simulation."""

    # TODO: making this frozen does not freeze all the submodel!
    model_config = ConfigDict(extra="forbid", frozen=True)

    _constituent_mapping: Dict[ConstituentString, ConstituentString] = {}

    # The following fields correspond to the toplevel fields of the json file descirbing a simulation
    atom1: ModelAtom
    atom2: Optional[Union[ModelAtom, ConstituentString]] = None
    classical_light1: Optional[ModelClassicalLight] = None
    classical_light2: Optional[Union[ModelClassicalLight, ConstituentString]] = None

    interactions: Optional[ModelInteractions] = None
    numerics: ModelNumerics = Field(default_factory=ModelNumerics)
    overlaps: Optional[ModelOverlaps] = None

    @cached_property
    def constituents(self) -> Dict[ConstituentString, BaseModelConstituent]:
        """Return a dictionary of all not None constituents."""
        constituents = {
            "atom1": self.atom1,
            "atom2": self.atom2,
            "classical_light1": self.classical_light1,
            "classical_light2": self.classical_light2,
        }
        return {k: v for k, v in constituents.items() if v is not None}

    @cached_property
    def parameter_range_options(self) -> ParameterRangeOptions:
        """Return a collection of all parameter ranges."""
        parameters = {}
        for submodel in [self.interactions, self.atom1, self.atom2, self.classical_light1, self.classical_light2]:
            if submodel is None:
                continue
            for k, v in iter(submodel):
                if isinstance(v, ParameterRange):
                    parameters[k] = v
        return ParameterRangeOptions(parameters=parameters)

    @field_validator("atom2", "classical_light2", mode="after")
    @classmethod
    def validate_constituents(
        cls, const: Union[BaseModelConstituent, ConstituentString], info: ValidationInfo
    ) -> BaseModelConstituent:
        """Validate the constituents, i.e. replace constituents given as reference
        to another constituent via a string as reference to the same object.

        To keep track of which constituents are pointing to the same object
        we store the mappings in the dictionary _constituent_mapping.
        This can then be used to check, wether we use the same atoms/basis,
        as well as to also dump the model again as reference via a string (see self.).
        """
        if const is None or isinstance(const, BaseModelConstituent):
            return const

        # else const is a string/reference to another constituent
        const_name = info.field_name
        new_const = info.data.get(const, None)
        if new_const is None:
            raise ValueError(f"{const_name} given as reference to {const}, but {const} is not defined.")
        info.data.setdefault("_constituent_mapping", {})[const_name] = const
        return new_const

    @field_validator("interactions", "overlaps", mode="before")
    @classmethod
    def validate_combined_states(cls, submodel: Optional[dict], info: ValidationInfo) -> dict:
        """Validate combined_states, i.e. replace states, that are given as reference
        to a state of a constituent.state_of_interest (via an index) with the actual state.
        """
        if submodel is None:
            return None

        combined_states = submodel.get("combined_states_of_interest", [])
        for csoi in combined_states:
            for constit, state in csoi.items():
                constituent = info.data.get(constit, None)
                if constituent is None:
                    raise ValueError(f"Key {constit} in combined_states_of_interest is not a valid constituent.")
                if not isinstance(state, int):  # state is already a BaseModelState
                    continue
                # else: state is an index
                states_of_interest = constituent.states_of_interest
                if state > len(states_of_interest) or state < 0:
                    raise ValueError(
                        f"Index {state} in combined for {constit} is out of range of {constit}.states_of_interest."
                    )
                csoi[constit] = states_of_interest[state]
        return submodel

    @field_serializer("atom2", "classical_light2", mode="wrap")
    def serialize_constituents(
        self, constituent: BaseModelConstituent, nxt: SerializerFunctionWrapHandler, info: FieldSerializationInfo
    ) -> Union[BaseModelConstituent, ConstituentString]:
        """Serialize the constituents, i.e. if the constituent is a reference to another constituent,
        return the name of the other constituent (this is done using self._constituent_mapping).
        """
        name = info.field_name
        if name in self._constituent_mapping:
            return self._constituent_mapping[name]
        return nxt(constituent, info)

    # Sanity checks
    @model_validator(mode="after")
    def sanity_check_fields(self) -> "ModelSimulation":
        """Check wether all atoms have the same applied fields, if not raise a warning."""
        if self.atom2 is None:
            return self

        for field in ["efield_x", "efield_y", "efield_z", "bfield_x", "bfield_y", "bfield_z"]:
            if getattr(self.atom1, field) != getattr(self.atom2, field):
                # TODO use logger or warnings
                print(f"Warning: {field} is different for atom1 and atom2, are you sure this is what you want?")

        return self

    @model_validator(mode="after")
    def compute_properties(self) -> "ModelSimulation":
        """Compute properties of the model, so also there it is checked if everything works without errors.

        This is especially important for parameter_range_options.
        """
        self.constituents  # noqa: B018
        self.parameter_range_options  # noqa: B018
        return self

    @classmethod
    def model_validate_json_file(cls, path: str, **kwargs) -> "ModelSimulation":
        """Validate the model from a json file."""
        with open(path, encoding="utf-8") as f:
            return cls.model_validate_json(f.read(), **kwargs)

    def model_dump_json_file(self, path: str, **kwargs) -> None:
        """Dump the model to a json file."""
        with open(path, "w", encoding="utf-8") as f:
            f.write(self.model_dump_json(**kwargs))
