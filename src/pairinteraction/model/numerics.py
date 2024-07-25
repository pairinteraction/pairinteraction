"""Pydantic models for numerical settings."""

from pathlib import Path
from typing import Literal, Optional, Union

from pydantic import BaseModel, ConfigDict, Field, field_serializer, field_validator

from pairinteraction.model.types.simple_types import Vector

QuantizationAxisString = Literal["x", "y", "z", "efield", "bfield", "distance"]


class ModelNumerics(BaseModel):
    """Model for various numerical parameters."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    quantization_axis: Union[QuantizationAxisString, Vector] = Field(default=[0, 0, 1])
    use_diamagnetism: bool = False
    rydberg_rydberg_multipole_order: int = Field(default=3, ge=3)

    database_directory: Optional[str] = None  # TODO should we check if it exists (using DirectoryPath) or not?
    precision_sparse_matrices: float = Field(default=1e-6, ge=0, le=1)
    fuzziness_in_database_search: float = Field(default=1e-1, ge=0)

    diagonalization_method: Literal[
        "lapacke", "feast_dense", "feast_sparse", "auto", "cuda"
    ] = "lapacke"  # TODO NEW: remove cuda for new version
    diagonalization_eigenspace_size: Optional[int] = Field(default=None, ge=1)

    path_cache: Optional[str] = None  # TODO NEW: remove for new version
    path_quantum_defects_database: Optional[str] = None  # TODO NEW: remove for new version
    diagonalize_threshold: float = 1e-10  # TODO NEW: remove for new version
    radial_wavefunction_method: Literal["numerov", "whittaker"] = "numerov"  # TODO NEW: remove for new version

    @field_validator("quantization_axis", mode="before")
    @classmethod
    def validate_quantization_axis(cls, axis: Union[QuantizationAxisString, Vector]) -> Vector:
        """Validate the quantization axis, i.e. if it is given as string (referencing some direction)
        convert it to a list, representing the direction.
        """
        if isinstance(axis, list):
            return axis

        if axis in ["x", "y", "z"]:
            axis = {"x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1]}[axis]
            return axis

        if axis in ["efield", "bfield"]:
            raise NotImplementedError(f"Quantization axis '{axis}' not yet implemented.")
            # TODO how to get the efield and bfield in numerics, or put this as model_validator for ModelSimulation?

        raise ValueError(f"Quantization axis '{axis}' not understood.")

    @field_serializer("database_directory", "path_cache")
    def serialize_path(self, path: Optional[Path]) -> Optional[str]:
        """Serialize paths as strings."""
        if path is None:
            return None
        return str(path)
