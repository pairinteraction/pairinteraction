"""Class for handling parameters."""

import types
from abc import ABC, abstractmethod
from typing import Any, Dict, Generic, List, Optional, Tuple, Type, TypeVar, Union, get_args

import numpy as np
from pydantic import GetCoreSchemaHandler
from pydantic_core import core_schema
from typing_extensions import Self

from pairinteraction.model.types.simple_types import Symmetry

ParamType = TypeVar("ParamType")
RawType = TypeVar("RawType")


class BaseParameter(ABC, Generic[ParamType, RawType]):
    """BaseParameter class."""

    def __init__(self, raw: RawType):
        self.raw = raw

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.raw})"

    @classmethod
    def __get_pydantic_core_schema__(cls, source: Type[Any], handler: GetCoreSchemaHandler) -> core_schema.CoreSchema:
        return core_schema.no_info_after_validator_function(
            cls._validate_raw,
            cls._get_pydantic_raw_schema(),
            serialization=core_schema.plain_serializer_function_ser_schema(
                cls._serialize,
                info_arg=False,
                return_schema=cls._get_pydantic_raw_schema(),
            ),
        )

    @classmethod
    @abstractmethod
    def _get_pydantic_raw_schema(cls) -> core_schema.CoreSchema:
        """Get the pydantic schema for the raw data."""

    @classmethod
    @abstractmethod
    def _get_pydantic_parameter_schema(cls) -> core_schema.CoreSchema:
        """Get the pydantic schema for the type of a single parameter value."""

    @classmethod
    def _validate_raw(cls, raw: RawType) -> Self:
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
    def get_list(self) -> Optional[List[ParamType]]:
        """Get the list of values (None for a constant parameter)."""

    @abstractmethod
    def get_min(self) -> ParamType:
        """Get the minimum value."""

    @abstractmethod
    def get_max(self) -> ParamType:
        """Get the maximum value."""


class BaseParameterIterable(BaseParameter[ParamType, RawType]):
    """Baseclass for an iterable parameter (e.g. list or range)."""

    list = None

    def __init__(self, raw: RawType):
        assert getattr(self, "list", None) is not None, "BaseParameterIterable needs to have a list attribute"
        super().__init__(raw)

    def get_value(self, step: int) -> ParamType:
        return self.list[step]

    def get_size(self) -> int:
        return len(self.list)

    def get_list(self) -> List[ParamType]:
        return self.list

    def get_min(self) -> ParamType:
        return min(self.list)

    def get_max(self) -> ParamType:
        return max(self.list)


class ParameterConstant(BaseParameter[ParamType, ParamType]):
    """Class for a constant parameter."""

    @classmethod
    def _get_pydantic_raw_schema(cls) -> core_schema.CoreSchema:
        return cls._get_pydantic_parameter_schema()

    def __init__(self, value: ParamType):
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


class ParameterList(BaseParameterIterable[ParamType, List[ParamType]]):
    """Class for a list of parameters."""

    @classmethod
    def _get_pydantic_raw_schema(cls) -> core_schema.CoreSchema:
        return core_schema.list_schema(cls._get_pydantic_parameter_schema(), min_length=1)

    def __init__(self, values: Union[List[ParamType], Tuple[ParamType]]):
        self.list = values
        super().__init__(values)


class ParameterRange(BaseParameterIterable[ParamType, Dict[str, Union[ParamType, int]]]):
    """Class for a parameter range."""

    @classmethod
    def _get_pydantic_raw_schema(cls) -> core_schema.CoreSchema:
        # TODO enforce check that start and stop are of type ParamType
        return core_schema.dict_schema(keys_schema=core_schema.str_schema(), min_length=3, max_length=3)

    def __init__(self, start: ParamType, stop: ParamType, steps: int):
        self.list = list(np.linspace(start, stop, steps))
        self.raw = {"start": start, "stop": stop, "steps": steps}

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(" + ", ".join([f"{k}={v!r}" for k, v in self.raw.items()]) + ")"

    @classmethod
    def _validate_raw(cls, raw: Dict[str, Union[ParamType, int]]) -> Self:
        if not all(key in raw for key in ["start", "stop", "steps"]):
            raise ValueError("ParameterRange needs to have keys 'start', 'stop', 'steps'")
        return cls(raw["start"], raw["stop"], raw["steps"])


class ParameterInt(BaseParameter):
    @classmethod
    def _get_pydantic_parameter_schema(cls) -> core_schema.CoreSchema:
        return core_schema.int_schema()


class ParameterFloat(BaseParameter):
    @classmethod
    def _get_pydantic_parameter_schema(cls) -> core_schema.CoreSchema:
        return core_schema.float_schema()


class ParameterSymmetry(BaseParameter):
    @classmethod
    def _get_pydantic_parameter_schema(cls) -> core_schema.CoreSchema:
        return core_schema.literal_schema(get_args(Symmetry))


ParameterConstantInt = types.new_class("ParameterConstantInt", (ParameterConstant[int], ParameterInt), {})
ParameterListInt = types.new_class("ParameterListInt", (ParameterList[int], ParameterInt), {})
ParameterRangeInt = types.new_class("ParameterRangeInt", (ParameterRange[int], ParameterInt), {})

ParameterConstantFloat = types.new_class("ParameterConstantFloat", (ParameterConstant[float], ParameterFloat), {})
ParameterListFloat = types.new_class("ParameterListFloat", (ParameterList[float], ParameterFloat), {})
ParameterRangeFloat = types.new_class("ParameterRangeFloat", (ParameterRange[float], ParameterFloat), {})

ParameterConstantSymmetry = types.new_class(
    "ParameterConstantSymmetry", (ParameterConstant[Symmetry], ParameterSymmetry), {}
)
ParameterListSymmetry = types.new_class("ParameterListSymmetry", (ParameterList[Symmetry], ParameterSymmetry), {})
