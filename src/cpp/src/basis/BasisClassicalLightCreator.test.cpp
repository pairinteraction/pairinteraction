#include "pairinteraction/basis/BasisClassicalLightCreator.hpp"

#include "pairinteraction/basis/BasisClassicalLight.hpp"

#include <doctest/doctest.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("create a classical light basis") {
    float energy{0.03};
    auto basis = BasisClassicalLightCreator<float>()
                     .set_photon_energy(energy)
                     .restrict_quantum_number_q(-3, 3)
                     .create();
    for (const auto &ket : *basis) {
        DOCTEST_CHECK(ket->get_photon_energy() == energy);
        DOCTEST_CHECK(ket->get_energy() == ket->get_photon_energy() * ket->get_quantum_number_q());
        DOCTEST_MESSAGE("Ket: ", *ket);
    }
}
} // namespace pairinteraction
