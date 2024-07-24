#include "pintr/ket/KetClassicalLightCreator.hpp"

#include "pintr/ket/KetClassicalLight.hpp"
#include "pintr/utils/streamed.hpp"

#include <doctest/doctest.h>
#include <spdlog/spdlog.h>

namespace pintr {
DOCTEST_TEST_CASE("create a classical light ket") {

    float photon_energy{0.02};
    int q{1};
    auto lightket = KetClassicalLightCreator<float>(photon_energy, q).create();

    DOCTEST_CHECK(lightket->get_quantum_number_q() == q);
    DOCTEST_CHECK(lightket->get_photon_energy() == photon_energy);
    DOCTEST_CHECK(lightket->get_energy() == q * photon_energy);
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Ket: {}", fmt::streamed(*lightket));
}
} // namespace pintr
