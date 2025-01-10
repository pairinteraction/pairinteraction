"""Class for handling parameters."""

import typing
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, Generic, Literal, Optional, TypeVar, Union

import numpy as np
from pydantic import GetCoreSchemaHandler
from pydantic_core import core_schema

if TYPE_CHECKING:
    from typing_extensions import Self

ParamType = TypeVar("ParamType")
RawType = TypeVar("RawType")


class BaseParameter(ABC, Generic[ParamType, RawType]):
    """BaseParameter class."""

    def __init__(self, raw: RawType) -> None:
        self.raw = raw

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.raw})"

    @classmethod
    def __get_pydantic_core_schema__(cls, source: type[Any], handler: GetCoreSchemaHandler) -> core_schema.CoreSchema:
        return core_schema.no_info_after_validator_function(
            cls._validate_raw,
            cls._get_pydantic_raw_schema(source),
            serialization=core_schema.plain_serializer_function_ser_schema(
                cls._serialize,
                info_arg=False,
                return_schema=cls._get_pydantic_raw_schema(source),
            ),
        )

    @classmethod
    @abstractmethod
    def _get_pydantic_raw_schema(cls, source: type[Any]) -> core_schema.CoreSchema:
        """Get the pydantic schema for the raw data."""

    @classmethod
    def _get_pydantic_parameter_schema(cls, source: type[Any]) -> core_schema.CoreSchema:
        """Get the pydantic schema for the type of a single parameter value."""
        generic_type = typing.get_args(source)
        if len(generic_type) == 0:
            raise ValueError(f"Unsupported type {source}")
        generic_type = generic_type[0]
        if isinstance(generic_type, type):
            if issubclass(generic_type, int):
                return core_schema.int_schema()
            if issubclass(generic_type, float):
                return core_schema.float_schema()
        if typing.get_origin(generic_type) == Literal:
            return core_schema.literal_schema(typing.get_args(generic_type))
        raise ValueError(f"Unsupported type {source}")

    @classmethod
    def _validate_raw(cls, raw: RawType) -> "Self":
        """Validate the raw input."""
        return cls(raw)

    @staticmethod
    def _serialize(parameter: "BaseParameter") -> RawType:
        return parameter.dump()

    def dump(self) -> RawType:
        """Dump the raw parameter."""
        return self.raw

    @abstractmethod
    def get_value(self, step: int) -> ParamType:
        """Get the value at a specific step."""

    @abstractmethod
    def get_size(self) -> Optional[int]:
        """Get the size of the parameter (None for a constant parameter)."""

    @abstractmethod
    def get_list(self) -> Optional[list[ParamType]]:
        """Get the list of values (None for a constant parameter)."""

    @abstractmethod
    def get_min(self) -> ParamType:
        """Get the minimum value."""

    @abstractmethod
    def get_max(self) -> ParamType:
        """Get the maximum value."""


class BaseParameterIterable(BaseParameter[ParamType, RawType]):
    """Baseclass for an iterable parameter (e.g. list or range)."""

    _list = None

    def __init__(self, raw: RawType) -> None:
        assert getattr(self, "_list", None) is not None, "BaseParameterIterable needs to have a _list attribute"
        super().__init__(raw)

    def get_value(self, step: int) -> ParamType:
        return self._list[step]

    def get_size(self) -> int:
        return len(self._list)

    def get_list(self) -> list[ParamType]:
        return self._list

    def get_min(self) -> ParamType:
        return min(self._list)

    def get_max(self) -> ParamType:
        return max(self._list)


class ParameterConstant(BaseParameter[ParamType, ParamType]):
    """Class for a constant parameter."""

    @classmethod
    def _get_pydantic_raw_schema(cls, source: type[Any]) -> core_schema.CoreSchema:
        return cls._get_pydantic_parameter_schema(source)

    def __init__(self, value: ParamType) -> None:
        self.value = value
        super().__init__(value)

    def get_value(self, step: Optional[int] = None) -> ParamType:
        return self.value

    def get_size(self) -> None:
        return None

    def get_list(self) -> None:
        return None

    def get_min(self) -> ParamType:
        return self.value

    def get_max(self) -> ParamType:
        return self.value


class ParameterList(BaseParameterIterable[ParamType, list[ParamType]]):
    """Class for a list of parameters."""

    @classmethod
    def _get_pydantic_raw_schema(cls, source: type[Any]) -> core_schema.CoreSchema:
        return core_schema.list_schema(cls._get_pydantic_parameter_schema(source), min_length=1)

    def __init__(self, values: Union[list[ParamType], tuple[ParamType]]) -> None:
        self._list = values
        super().__init__(values)


class ParameterRange(BaseParameterIterable[ParamType, dict[str, Union[ParamType, int]]]):
    """Class for a parameter range."""

    @classmethod
    def _get_pydantic_raw_schema(cls, source: type[Any]) -> core_schema.CoreSchema:
        # TODO enforce check that start and stop are of type ParamType
        return core_schema.dict_schema(keys_schema=core_schema.str_schema(), min_length=3, max_length=3)

    def __init__(self, start: ParamType, stop: ParamType, steps: int) -> None:
        self._list = list(np.linspace(start, stop, steps))
        self.raw = {"start": start, "stop": stop, "steps": steps}

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(" + ", ".join([f"{k}={v!r}" for k, v in self.raw.items()]) + ")"

    @classmethod
    def _validate_raw(cls, raw: dict[str, Union[ParamType, int]]) -> "Self":
        if not all(key in raw for key in ["start", "stop", "steps"]):
            raise ValueError("ParameterRange needs to have keys 'start', 'stop', 'steps'")
        return cls(raw["start"], raw["stop"], raw["steps"])
