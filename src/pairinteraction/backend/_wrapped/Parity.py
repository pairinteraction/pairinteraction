from typing import Literal, get_args

import pairinteraction.backend._backend as _backend

Parity = Literal["EVEN", "ODD", "UNKNOWN"]

CPPParity = _backend.Parity


def get_cpp_parity(parity: Parity) -> CPPParity:
    if parity.upper() not in get_args(Parity):
        raise ValueError(f"Unknown parity '{parity}', should be one of {Parity}")
    try:
        return getattr(CPPParity, f"{parity.upper()}")()
    except AttributeError as err:
        raise ValueError(f"Unknown parity 'Parity.{parity.upper()}'") from err
