from pairinteraction_backend import KetClassicalLightCreatorFloat

def test_ket_classical_light_creator():
    ket = KetClassicalLightCreatorFloat(1.0, 2).create()
    assert ket.get_photon_energy() == 1.0
    assert ket.get_quantum_number_q() == 2
