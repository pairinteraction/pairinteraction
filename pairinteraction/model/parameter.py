"""Class for handling parameters."""

from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Tuple, Type, Union

import numpy as np
from pydantic import GetCoreSchemaHandler
from pydantic_core import core_schema


class BaseParameter(ABC):
    """BaseParameter class."""

    _pydantic_schema = core_schema.any_schema()

    def __init__(self, raw: Any):
        self.raw = raw

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.raw})"

    @classmethod
    def __get_pydantic_core_schema__(cls, source: Type[Any], handler: GetCoreSchemaHandler) -> core_schema.CoreSchema:
        return core_schema.no_info_after_validator_function(
            cls._validate_raw,
            cls._pydantic_schema,
            serialization=core_schema.plain_serializer_function_ser_schema(
                cls._serialize,
                info_arg=False,
                return_schema=cls._pydantic_schema,
            ),
        )

    @classmethod
    def _validate_raw(cls, raw: Any) -> Union["ParameterConstant", "ParameterList", "ParameterRange"]:
        """Validate the raw input."""
        return cls(raw)

    @staticmethod
    def _serialize(parameter: "BaseParameter") -> Any:
        return parameter.dump()

    def dump(self):
        """Dump the raw parameter."""
        return self.raw

    @abstractmethod
    def get_value(self, step: int) -> float:
        """Get the value at a specific step."""

    @abstractmethod
    def get_size(self) -> Optional[int]:
        """Get the size of the parameter (None for a constant parameter)."""

    @abstractmethod
    def get_list(self) -> Optional[List[float]]:
        """Get the list of values (None for a constant parameter)."""

    @abstractmethod
    def get_min(self) -> float:
        """Get the minimum value."""

    @abstractmethod
    def get_max(self) -> float:
        """Get the maximum value."""


class ParameterConstant(BaseParameter):
    """Class for a constant parameter."""

    _pydantic_schema = core_schema.float_schema()

    def __init__(self, value: float):
        self.value = value
        super().__init__(value)

    def get_value(self, step: Optional[int] = None) -> float:
        return self.value

    def get_size(self) -> None:
        return None

    def get_list(self) -> None:
        return None

    def get_min(self) -> float:
        return self.value

    def get_max(self) -> float:
        return self.value


class ParameterList(BaseParameter):
    """Class for a list of parameters."""

    _pydantic_schema = core_schema.list_schema(core_schema.float_schema(), min_length=1)

    def __init__(self, values: Union[List[float], Tuple[float]]):
        self.list = values
        super().__init__(values)

    def get_value(self, step: int) -> float:
        return self.list[step]

    def get_size(self) -> int:
        return len(self.list)

    def get_list(self) -> List[float]:
        return self.list

    def get_min(self) -> float:
        return min(self.list)

    def get_max(self) -> float:
        return max(self.list)


class ParameterRange(ParameterList):
    """Class for a parameter range."""

    _pydantic_schema = core_schema.dict_schema(keys_schema=core_schema.str_schema())

    def __init__(self, start: float, stop: float, steps: int):
        self.list = list(np.linspace(start, stop, steps))
        self.raw = {"start": start, "stop": stop, "steps": steps}

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(" + ", ".join([f"{k}={v!r}" for k, v in self.raw.items()]) + ")"

    @classmethod
    def _validate_raw(cls, raw: Dict[str, Union[float, int]]) -> "ParameterRange":
        return cls(raw["start"], raw["stop"], raw["steps"])


UnionParameter = Union[ParameterConstant, ParameterList, ParameterRange]
