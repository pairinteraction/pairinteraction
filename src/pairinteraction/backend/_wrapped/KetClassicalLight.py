from typing import Optional, Union

import pairinteraction.backend._backend as _backend


class KetClassicalLightBase:
    _KetClassicalLightCreator: Union[
        type[_backend.KetClassicalLightCreatorFloat], type[_backend.KetClassicalLightCreatorDouble]
    ]

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

    @classmethod
    def _from_cpp_object(cls, cpp_obj: Union[_backend.KetClassicalLightFloat, _backend.KetClassicalLightDouble]):
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    # # # Define convenience properties # # #
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
