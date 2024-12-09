#include "pairinteraction/ket/KetClassicalLightCreator.hpp"

#include "pairinteraction/ket/KetClassicalLight.hpp"

#include <doctest/doctest.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("create a classical light ket") {

    float photon_energy{0.02};
    int q{1};
    auto lightket = KetClassicalLightCreator<float>(photon_energy, q).create();

    DOCTEST_CHECK(lightket->get_quantum_number_q() == q);
    DOCTEST_CHECK(lightket->get_photon_energy() == photon_energy);
    DOCTEST_CHECK(lightket->get_energy() == q * photon_energy);
    DOCTEST_MESSAGE("Ket: ", *lightket);
}
} // namespace pairinteraction
