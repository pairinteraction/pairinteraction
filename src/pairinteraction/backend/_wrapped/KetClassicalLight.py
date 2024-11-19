from typing import Optional, Union

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.Ket import KetBase

UnionCPPKetClassicalLight = Union[_backend.KetClassicalLightFloat, _backend.KetClassicalLightDouble]
UnionTypeCPPKetAtomCreator = Union[
    type[_backend.KetClassicalLightCreatorFloat], type[_backend.KetClassicalLightCreatorDouble]
]


class KetClassicalLightBase(KetBase):
    _cpp: UnionCPPKetClassicalLight
    _KetClassicalLightCreator: UnionTypeCPPKetAtomCreator

    def __init__(
        self,
        photon_energy: Optional[float] = None,
        q: Optional[int] = None,
    ) -> None:
        creator = self._KetClassicalLightCreator()
        if photon_energy is not None:
            creator.set_photon_energy(photon_energy)
        if q is not None:
            creator.set_quantum_number_q(q)
        self._cpp = creator.create()

    @property
    def q(self) -> int:
        return self._cpp.get_quantum_number_q()

    @property
    def photon_energy(self) -> float:
        return self._cpp.get_photon_energy()


class KetClassicalLightFloat(KetClassicalLightBase):
    _cpp: _backend.KetClassicalLightFloat
    _KetClassicalLightCreator = _backend.KetClassicalLightCreatorFloat


class KetClassicalLightDouble(KetClassicalLightBase):
    _cpp: _backend.KetClassicalLightDouble
    _KetClassicalLightCreator = _backend.KetClassicalLightCreatorDouble
