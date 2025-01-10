"""Pydantic model for the complete simulation."""

from functools import cached_property
from typing import TYPE_CHECKING, Optional, Union

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    FieldSerializationInfo,
    SerializerFunctionWrapHandler,
    field_serializer,
    model_validator,
)

if TYPE_CHECKING:
    from typing_extensions import Self

from pairinteraction.model.constituents.atom import ModelAtomMQDT, ModelAtomSQDT
from pairinteraction.model.constituents.base import BaseModelConstituent
from pairinteraction.model.constituents.classical_light import ModelClassicalLight
from pairinteraction.model.interactions import ModelInteractions
from pairinteraction.model.numerics import ModelNumerics
from pairinteraction.model.overlaps import ModelOverlaps
from pairinteraction.model.states.base import BaseModelState
from pairinteraction.model.types.parameter import BaseParameterIterable
from pairinteraction.model.types.simple_types import ConstituentString, HalfInt

UnionModelAtom = Union[ModelAtomSQDT[int], ModelAtomSQDT[HalfInt], ModelAtomMQDT[int], ModelAtomMQDT[HalfInt]]


class ModelSimulation(BaseModel):
    """Model for the input of a pairinteraction simulation."""

    model_config = ConfigDict(extra="forbid", frozen=False, validate_assignment=True)

    # The following fields correspond to the toplevel fields of the json file descirbing a simulation
    atom1: UnionModelAtom
    atom2: Optional[Union[UnionModelAtom, ConstituentString]] = None
    classical_light1: Optional[ModelClassicalLight] = None
    classical_light2: Optional[Union[ModelClassicalLight, ConstituentString]] = None

    interactions: Optional[ModelInteractions] = None
    numerics: ModelNumerics = Field(default_factory=ModelNumerics)
    overlaps: Optional[ModelOverlaps] = None

    @property
    def constituents(self) -> dict[ConstituentString, BaseModelConstituent]:
        """Return a dictionary of all not None constituents."""
        constituents = {
            "atom1": self.atom1,
            "atom2": self.atom2,
            "classical_light1": self.classical_light1,
            "classical_light2": self.classical_light2,
        }
        return {k: v for k, v in constituents.items() if v is not None}

    @property
    def unique_constituents(self) -> dict[ConstituentString, BaseModelConstituent]:
        """Return a dictionary of all unique constituents."""
        constituents = {
            "atom1": self.atom1,
            "atom2": self.atom2,
            "classical_light1": self.classical_light1,
            "classical_light2": self.classical_light2,
        }
        return {k: v for k, v in constituents.items() if v is not None and k not in self._constituent_mapping}

    @cached_property
    def dict_of_parameter_lists(self) -> dict[str, BaseParameterIterable]:
        """Return a collection of all parameter ranges."""
        parameters = {}
        for submodel in [self.interactions, self.atom1, self.atom2, self.classical_light1, self.classical_light2]:
            if submodel is None:
                continue
            for k, v in iter(submodel):
                if isinstance(v, BaseParameterIterable):
                    parameters[k] = v
        return parameters

    @cached_property
    def parameter_size(self) -> int:
        """Return the number of steps of the parameter lists, and make sure they are all the same size."""
        all_sizes = [p.get_size() for p in self.dict_of_parameter_lists.values()]
        if len(all_sizes) == 0:
            return 0
        if not all(s == all_sizes[0] for s in all_sizes):
            raise ValueError("All parameter lists must have the same size.")
        return all_sizes[0]

    @model_validator(mode="after")
    def validate_constituents_references(self) -> "Self":
        """Validate the constituents.

        This means replace constituents given as reference
        to another constituent via a string as reference to the same object.

        To keep track of which constituents are pointing to the same object
        we store the mappings in the dictionary _constituent_mapping.
        This can then be used to check, wether we use the same atoms/basis,
        as well as to also dump the model again as reference via a string (see self.serialize_constituents).
        """
        if not hasattr(self, "_constituent_mapping"):
            self._constituent_mapping = {}
        for const_name in ["atom", "classical_light"]:
            if isinstance(getattr(self, const_name + "2"), str):
                self._constituent_mapping[const_name + "2"] = const_name + "1"
                setattr(self, const_name + "2", getattr(self, const_name + "1"))
        return self

    @model_validator(mode="after")
    def validate_submodel_combined_states(self) -> "Self":
        """Validate combined_states of the submodel 'interactions' and 'overlaps'.

        This means replace states, that are given as reference
        to a state of a constituent.state_of_interest (via an index) with the actual state.
        """
        for submodel_name in ["interactions", "overlaps"]:
            submodel = getattr(self, submodel_name)
            if submodel is None:
                continue
            for csoi in submodel.combined_states_of_interest:
                for constit, state in csoi.items():
                    constituent = getattr(self, constit)
                    if constituent is None:
                        raise ValueError(
                            f"{submodel_name}: Key {constit} in combined_states_of_interest is not a valid constituent."
                        )
                    if isinstance(state, BaseModelState):
                        continue
                    # else: state is an index
                    if state > len(constituent.states_of_interest) or state < 0:
                        raise ValueError(
                            f"Index {state} in {submodel_name}.combined_states_of_interest "
                            f"is out of range of {constit}.states_of_interest."
                        )
                    csoi[constit] = constituent.states_of_interest[state]

        return self

    @model_validator(mode="after")
    def sanity_check_fields(self) -> "Self":
        """Check wether all atoms have the same applied fields, if not raise a warning."""
        if self.atom2 is None:
            return self

        for field in ["efield_x", "efield_y", "efield_z", "bfield_x", "bfield_y", "bfield_z"]:
            if getattr(self.atom1, field) != getattr(self.atom2, field):
                # TODO use logger or warnings
                print(f"Warning: {field} is different for atom1 and atom2, are you sure this is what you want?")

        return self

    @model_validator(mode="after")
    def compute_properties(self) -> "Self":
        """Compute properties of the model, so also there it is checked if everything works without errors.

        This is especially important for dict_of_parameter_lists.
        """
        self.constituents  # noqa: B018
        self.dict_of_parameter_lists  # noqa: B018
        self.parameter_size  # noqa: B018
        return self

    @field_serializer("atom2", "classical_light2", mode="wrap")
    def serialize_constituents(
        self,
        constituent: Optional[BaseModelConstituent],
        nxt: SerializerFunctionWrapHandler,
        info: FieldSerializationInfo,
    ) -> Optional[Union[dict, ConstituentString]]:
        """Serialize the constituents.

        If the constituent is a reference to another constituent
        return the name of the other constituent (this is done using self._constituent_mapping).
        """
        name = info.field_name
        if name in self._constituent_mapping:
            return self._constituent_mapping[name]
        return nxt(constituent, info)

    @classmethod
    def model_validate_json_file(cls, path: str, **kwargs) -> "Self":
        """Validate the model from a json file."""
        with open(path, encoding="utf-8") as f:
            return cls.model_validate_json(f.read(), **kwargs)

    def model_dump_json_file(self, path: str, **kwargs) -> None:
        """Dump the model to a json file."""
        with open(path, "w", encoding="utf-8") as f:
            f.write(self.model_dump_json(**kwargs))
