from typing import Literal, Union

from pairinteraction.backend import _backend

Diagonalizer = Literal["Eigen", "Lapacke", "Feast"]

UnionCPPDiagonalizer = Union[
    _backend.DiagonalizerInterfaceFloat,
    _backend.DiagonalizerInterfaceComplexFloat,
    _backend.DiagonalizerInterfaceDouble,
    _backend.DiagonalizerInterfaceComplexDouble,
]
