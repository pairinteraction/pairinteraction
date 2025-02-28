import pairinteraction.real as pi


def test_ket() -> None:
    ket = pi.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    assert ket.species == "Rb"
    assert ket.n == 60
    assert ket.l == 0
    assert ket.j == 0.5
    assert ket.m == 0.5
