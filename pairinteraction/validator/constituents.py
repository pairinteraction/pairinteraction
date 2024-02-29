"""Pydantic models for constituents."""

from typing import List, Optional, Union

import numpy as np
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    ValidationInfo,
    field_validator,
    model_validator,
)

from pairinteraction.validator import misc
from pairinteraction.validator.misc import ExtraField, HalfInt, LeftToRightField, SpeciesString
from pairinteraction.validator.parameter import UnionParameter
from pairinteraction.validator.reusing_validators import (
    one_use_delta_and_soi,
    use_parameter_if_float,
)
from pairinteraction.validator.states import UnionModelStates


class BaseModelConstituent(BaseModel):
    """Pydantic base model corresponding to a System... ."""

    # TODO: is there a option to set frozen to True after all model_validators have been run?
    model_config = ConfigDict(extra="forbid", frozen=False, validate_assignment=True)


class ModelAtom(BaseModelConstituent):
    """Pydantic model corresponding to SystemAtom."""

    species: SpeciesString

    min_n: Union[int, float, None] = LeftToRightField(None)
    max_n: Union[int, float, None] = LeftToRightField(None)
    min_nu: Optional[float] = None
    max_nu: Optional[float] = None
    min_l: Union[int, float, None] = LeftToRightField(None)
    max_l: Union[int, float, None] = LeftToRightField(None)
    min_s: Union[HalfInt, float, None] = LeftToRightField(None)
    max_s: Union[HalfInt, float, None] = LeftToRightField(None)
    min_j: Union[HalfInt, float, None] = LeftToRightField(None)
    max_j: Union[HalfInt, float, None] = LeftToRightField(None)
    min_f: Optional[HalfInt] = None
    max_f: Optional[HalfInt] = None
    min_m: Optional[HalfInt] = None
    max_m: Optional[HalfInt] = None

    min_energy: Optional[float] = None
    max_energy: Optional[float] = None
    min_energy_after_diagonalization: Optional[float] = None
    max_energy_after_diagonalization: Optional[float] = None

    efield_x: UnionParameter = Field(0, validate_default=True)
    efield_y: UnionParameter = Field(0, validate_default=True)
    efield_z: UnionParameter = Field(0, validate_default=True)
    bfield_x: UnionParameter = Field(0, validate_default=True)
    bfield_y: UnionParameter = Field(0, validate_default=True)
    bfield_z: UnionParameter = Field(0, validate_default=True)

    # Optional use delta_attr and states_of_interest to define min_/max_attr
    # dont use BaseModelState here, but specific UnionModelStates to enable automatic parsing
    _used_delta: bool = False
    states_of_interest: List[UnionModelStates] = []
    delta_n: Optional[int] = ExtraField()
    delta_nu: Optional[float] = ExtraField()
    delta_l: Optional[int] = ExtraField()
    delta_s: Optional[int] = ExtraField()
    delta_j: Optional[int] = ExtraField()
    delta_f: Optional[int] = ExtraField()
    delta_m: Optional[int] = ExtraField()
    delta_energy: Optional[float] = ExtraField()
    delta_energy_after_diagonalization: Optional[float] = ExtraField()

    additionally_included_states: List[UnionModelStates] = []

    use_parameter_if_float = field_validator(
        "efield_x", "efield_y", "efield_z", "bfield_x", "bfield_y", "bfield_z", mode="before"
    )(use_parameter_if_float)

    @property
    def is_real(self) -> bool:
        """Return True if the atom can be described with a real Hamiltonian (i.e. no efield or bfield in y direction).

        This defines wether to use pairinteraction.SystemOneReal or SystemOneComplex.
        """
        return not np.any(self.efield_y.list) and not np.any(self.bfield_y.list)

    @field_validator("species", mode="before")
    @classmethod
    def old_species_names(cls, species: SpeciesString) -> SpeciesString:
        """Translate the species to the old format."""
        # TODO NEW: remove this whole field validator
        return misc.NEW_TO_OLD_SPECIES.get(species, species)

    @field_validator("states_of_interest", "additionally_included_states", mode="before")
    @classmethod
    def set_states_species(cls, states: List, info: ValidationInfo) -> List:
        """If states_of_interest are provided, automatically set (and check)
        the species of all the states to the species provided to the ModelAtom.

        This allows to define the species only once and not for each state separately.
        """
        species = info.data.get("species")
        for soi in states:
            if not soi.setdefault("species", species) == species:
                raise ValueError("species of states must be the same as for the constituent")
        return states

    @model_validator(mode="after")
    def use_delta_and_soi(self) -> "ModelAtom":
        """If states_of_interest is provided together with delta_attr instead of min/max_attr
        convert the delta_attr to min/max_attr .
        """
        # FIXME: this is a ugly hack to avoid running use_delta again when setting
        # and thus revalidating the min/max values
        if self._used_delta:
            return self
        self._used_delta = True
        for attr in ["n", "nu", "l", "s", "j", "f", "m", "energy", "energy_after_diagonalization"]:
            one_use_delta_and_soi(self, attr)
        return self


class ModelClassicalLight(BaseModelConstituent):
    """Pydantic model corresponding to SystemClassicalLight."""

    # TODO defaults for frequency, pointing_vector, intensity_per_polarization?
    frequency: float
    pointing_vector: List[float]
    intensity_per_polarization: List[float]

    min_n_pi: Optional[int] = None
    max_n_pi: Optional[int] = None
    min_n_plus: Optional[int] = None
    max_n_plus: Optional[int] = None
    min_n_minus: Optional[int] = None
    max_n_minus: Optional[int] = None
    min_n_total: Optional[int] = None
    max_n_total: Optional[int] = None

    min_energy_after_diagonalization: Optional[float] = None
    max_energy_after_diagonalization: Optional[float] = None

    # dont use BaseModelState here, but specific UnionModelStates to enable automatic parsing
    states_of_interest: List[UnionModelStates] = []

    @field_validator("states_of_interest", mode="before")
    @classmethod
    def set_soi_attributes(cls, states_of_interest: List, info: ValidationInfo) -> List:
        """If states_of_interest are provided, automatically set (and check) the attributes
        (frequency, pointing_vector and intensity_per_polarization) of all the states.

        This allows to define the attributes only once and not for each state separately.
        """
        for soi in states_of_interest:
            for attr in ["frequency", "pointing_vector", "intensity_per_polarization"]:
                value = info.data.get(attr)
                assert (
                    soi.setdefault(attr, value) == value
                ), f"{attr} of states_of_interest must be the same as for the constituent"
        return states_of_interest


UnionModelConstituent = Union[ModelAtom, ModelClassicalLight]
