from typing import Union

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.SystemAtom import (
    SystemAtomBase,
    SystemAtomComplexDouble,
    SystemAtomComplexFloat,
    SystemAtomDouble,
    SystemAtomFloat,
)
from pairinteraction.backend._wrapped.SystemCombined import (
    SystemCombinedBase,
    SystemCombinedComplexDouble,
    SystemCombinedComplexFloat,
    SystemCombinedDouble,
    SystemCombinedFloat,
)


def diagonalize(
    systems: list[Union[SystemAtomBase, SystemCombinedBase]],
    diagonalizer: Union[_backend.DiagonalizerInterfaceFloat, _backend.DiagonalizerInterfaceDouble, None] = None,
    precision: int = 12,
    eigenvalue_range: Union[_backend.RangeFloat, _backend.RangeDouble, None] = None,
) -> list[Union[SystemAtomBase, SystemCombinedBase]]:
    cpp_systems = [s._cpp for s in systems]
    diagonalize_fct = get_diagonalize_fct(systems[0])

    if diagonalizer is None:
        diagonalizer = systems[0]._DefaultDiagonalizer()
    if eigenvalue_range is None:
        diagonalize_fct(cpp_systems, diagonalizer, precision)
    else:
        diagonalize_fct(cpp_systems, diagonalizer, precision, eigenvalue_range)
    for i, system in enumerate(systems):
        system._cpp = cpp_systems[i]
        system.update_basis()
    return systems


def get_diagonalize_fct(system):
    if isinstance(system, SystemAtomFloat):
        return _backend.diagonalizeSystemAtomFloat
    elif isinstance(system, SystemAtomComplexFloat):
        return _backend.diagonalizeSystemAtomComplexFloat
    elif isinstance(system, SystemAtomDouble):
        return _backend.diagonalizeSystemAtomDouble
    elif isinstance(system, SystemAtomComplexDouble):
        return _backend.diagonalizeSystemAtomComplexDouble
    elif isinstance(system, SystemCombinedFloat):
        return _backend.diagonalizeSystemCombinedFloat
    elif isinstance(system, SystemCombinedComplexFloat):
        return _backend.diagonalizeSystemCombinedComplexFloat
    elif isinstance(system, SystemCombinedDouble):
        return _backend.diagonalizeSystemCombinedDouble
    elif isinstance(system, SystemCombinedComplexDouble):
        return _backend.diagonalizeSystemCombinedComplexDouble
    else:
        raise ValueError(f"Unknown system type {type(system)}")
