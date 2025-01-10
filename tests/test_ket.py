import pairinteraction.backend.float as pi


def test_ket_classical_light_creator() -> None:
    ket = pi.KetClassicalLight(photon_energy=1.0, q=2)
    assert ket.photon_energy == 1.0
    assert ket.q == 2
