#include "pintr/basis/BasisClassicalLightCreator.hpp"

#include "pintr/basis/BasisClassicalLight.hpp"
#include "pintr/utils/streamed.hpp"

#include <doctest/doctest.h>
#include <spdlog/spdlog.h>

namespace pintr {
DOCTEST_TEST_CASE("create a classical light basis") {
    float energy{0.03};
    auto basis = BasisClassicalLightCreator<float>()
                     .set_photon_energy(energy)
                     .restrict_quantum_number_q(-3, 3)
                     .create();
    for (auto ket : *basis) {
        DOCTEST_CHECK(ket->get_photon_energy() == energy);
        DOCTEST_CHECK(ket->get_energy() == ket->get_photon_energy() * ket->get_quantum_number_q());
        SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Ket: {}", fmt::streamed(*ket));
    }
}
} // namespace pintr
