from typing import Literal

from pairinteraction.backend import _backend

Parity = Literal["EVEN", "ODD", "UNKNOWN"]

CPPParity = _backend.Parity
