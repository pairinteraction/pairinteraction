"""Class for handling parameters."""

from abc import ABC, abstractmethod
from typing import Any, Dict, Generic, List, Optional, Tuple, Type, TypeVar, Union, get_args

import numpy as np
from pydantic import GetCoreSchemaHandler
from pydantic_core import core_schema

from pairinteraction.model.types import Symmetry

T = TypeVar("T")


class BaseParameter(ABC, Generic[T]):
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
    def get_value(self, step: int) -> T:
        """Get the value at a specific step."""

    @abstractmethod
    def get_size(self) -> Optional[int]:
        """Get the size of the parameter (None for a constant parameter)."""

    @abstractmethod
    def get_list(self) -> Optional[List[T]]:
        """Get the list of values (None for a constant parameter)."""

    @abstractmethod
    def get_min(self) -> T:
        """Get the minimum value."""

    @abstractmethod
    def get_max(self) -> T:
        """Get the maximum value."""


class ParameterConstant(BaseParameter[T]):
    """Class for a constant parameter."""

    _pydantic_schema = core_schema.any_schema()  # Unfortunatly, we cannot use T here

    def __init__(self, value: T):
        self.value = value
        super().__init__(value)

    def get_value(self, step: Optional[int] = None) -> T:
        return self.value

    def get_size(self) -> None:
        return None

    def get_list(self) -> None:
        return None

    def get_min(self) -> T:
        return self.value

    def get_max(self) -> T:
        return self.value


class ParameterList(BaseParameter[T]):
    """Class for a list of parameters."""

    _pydantic_schema = core_schema.list_schema(core_schema.any_schema(), min_length=1)

    def __init__(self, values: Union[List[T], Tuple[T]]):
        self.list = values
        super().__init__(values)

    def get_value(self, step: int) -> T:
        return self.list[step]

    def get_size(self) -> int:
        return len(self.list)

    def get_list(self) -> List[T]:
        return self.list

    def get_min(self) -> T:
        return min(self.list)

    def get_max(self) -> T:
        return max(self.list)


class ParameterRange(ParameterList[T]):
    """Class for a parameter range."""

    _pydantic_schema = core_schema.dict_schema(keys_schema=core_schema.str_schema())

    def __init__(self, start: T, stop: T, steps: int):
        self.list = list(np.linspace(start, stop, steps))
        self.raw = {"start": start, "stop": stop, "steps": steps}

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(" + ", ".join([f"{k}={v!r}" for k, v in self.raw.items()]) + ")"

    @classmethod
    def _validate_raw(cls, raw: Dict[str, Union[T, int]]) -> "ParameterRange":
        return cls(raw["start"], raw["stop"], raw["steps"])


class ParameterConstantInt(ParameterConstant[int]):
    _pydantic_schema = core_schema.int_schema()


class ParameterConstantFloat(ParameterConstant[float]):
    _pydantic_schema = core_schema.float_schema()


class ParameterConstantSymmetry(ParameterConstant[Symmetry]):
    _pydantic_schema = core_schema.literal_schema(get_args(Symmetry))


class ParameterListInt(ParameterList[int]):
    _pydantic_schema = core_schema.list_schema(core_schema.int_schema())


class ParameterListFloat(ParameterList[float]):
    _pydantic_schema = core_schema.list_schema(core_schema.float_schema())


class ParameterListSymmetry(ParameterList[Symmetry]):
    _pydantic_schema = core_schema.list_schema(core_schema.literal_schema(get_args(Symmetry)))


class ParameterRangeInt(ParameterRange[int]):
    _pydantic_schema = core_schema.dict_schema(
        keys_schema=core_schema.str_schema(),
        values_schema=core_schema.int_schema(),
    )


class ParameterRangeFloat(ParameterRange[float]):
    _pydantic_schema = core_schema.dict_schema(
        keys_schema=core_schema.str_schema(),
    )


UnionParameterInt = Union[ParameterConstantInt, ParameterListInt, ParameterRangeInt]
UnionParameterFloat = Union[ParameterConstantFloat, ParameterListFloat, ParameterRangeFloat]
UnionParameterSymmetry = Union[ParameterConstantSymmetry, ParameterListSymmetry]
