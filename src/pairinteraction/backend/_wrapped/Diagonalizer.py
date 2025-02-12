from typing import Literal, Union

from pairinteraction.backend import _backend

Diagonalizer = Literal["Eigen", "Lapacke", "Feast"]

UnionCPPDiagonalizer = Union[_backend.DiagonalizerInterfaceReal, _backend.DiagonalizerInterfaceComplex]
