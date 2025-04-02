# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Any, Literal, Optional, Union

from pairinteraction import _backend, _wrapped
from pairinteraction._wrapped.enums import FloatType, get_cpp_float_type

if TYPE_CHECKING:
    from pairinteraction._wrapped.system.system import SystemBase


Diagonalizer = Literal["eigen", "lapacke_evd", "lapacke_evr", "feast"]
UnionCPPDiagonalizer = Union[_backend.DiagonalizerInterfaceReal, _backend.DiagonalizerInterfaceComplex]
UnionCPPDiagonalizerType = Union[type[_backend.DiagonalizerInterfaceReal], type[_backend.DiagonalizerInterfaceComplex]]

_DiagonalizerDict: dict[str, dict[Diagonalizer, UnionCPPDiagonalizerType]] = {
    "real": {
        "eigen": _backend.DiagonalizerEigenReal,
        "lapacke_evd": _backend.DiagonalizerLapackeEvdReal,
        "lapacke_evr": _backend.DiagonalizerLapackeEvrReal,
        "feast": _backend.DiagonalizerFeastReal,
    },
    "complex": {
        "eigen": _backend.DiagonalizerEigenComplex,
        "lapacke_evd": _backend.DiagonalizerLapackeEvdComplex,
        "lapacke_evr": _backend.DiagonalizerLapackeEvrComplex,
        "feast": _backend.DiagonalizerFeastComplex,
    },
}


def get_cpp_diagonalizer(
    diagonalizer: Diagonalizer,
    system: "SystemBase[Any]",
    float_type: FloatType,
    m0: Optional[int] = None,
) -> UnionCPPDiagonalizer:
    if diagonalizer == "feast" and m0 is None:
        raise ValueError("m0 must be specified for the 'feast' diagonalizer")
    if diagonalizer != "feast" and m0 is not None:
        raise ValueError("m0 must not be specified if the diagonalizer is not 'feast'")

    if isinstance(system, (_wrapped.SystemAtomReal, _wrapped.SystemPairReal)):
        type_ = "real"
    elif isinstance(system, (_wrapped.SystemAtomComplex, _wrapped.SystemPairComplex)):
        type_ = "complex"
    else:
        raise TypeError(
            f"system must be of type SystemAtomReal, SystemPairReal, SystemAtomComplex, or SystemPairComplex, "
            f"not {type(system)}"
        )

    try:
        diagonalizer_class = _DiagonalizerDict[type_][diagonalizer]
    except KeyError:
        raise ValueError(
            f"Unknown diagonalizer '{diagonalizer}', should be one of {list(_DiagonalizerDict.keys())}"
        ) from None

    cpp_float_type = get_cpp_float_type(float_type)
    if diagonalizer == "feast":
        return diagonalizer_class(m0=m0, float_type=cpp_float_type)  # type: ignore [call-arg]
    return diagonalizer_class(float_type=cpp_float_type)  # type: ignore [call-arg]
