#include "basis/BasisClassicalLightCreator.hpp"

#include "basis/BasisClassicalLight.hpp"
#include "ket/KetClassicalLightCreator.hpp"

template <typename Scalar>
BasisClassicalLightCreator<Scalar> &
BasisClassicalLightCreator<Scalar>::set_photon_energy(real_t value) {
    photon_energy.emplace(value);
    return *this;
}

template <typename Scalar>
BasisClassicalLightCreator<Scalar> &
BasisClassicalLightCreator<Scalar>::restrict_quantum_number_q(int min, int max) {
    min_quantum_number_q.emplace(min);
    max_quantum_number_q.emplace(max);
    return *this;
}

template <typename Scalar>
BasisClassicalLight<Scalar> BasisClassicalLightCreator<Scalar>::create() const {

    int extracted_min_q, extracted_max_q;
    real_t extracted_photon_energy;
    if (max_quantum_number_q.has_value() && min_quantum_number_q.has_value() &&
        photon_energy.has_value()) {
        extracted_min_q = min_quantum_number_q.value();
        extracted_max_q = max_quantum_number_q.value();
        extracted_photon_energy = photon_energy.value();
    } else if (!photon_energy.has_value()) {
        throw std::runtime_error("photon energy not specified!");
    } else {
        throw std::runtime_error("photon number not restricted!");
    }

    std::vector<std::shared_ptr<const ket_t>> kets;
    kets.reserve(extracted_max_q - extracted_min_q + 1);

    auto ket_creator =
        KetClassicalLightCreator<real_t>().set_photon_energy(extracted_photon_energy);

    // for (const int q : std::views::iota(extracted_min_q, extracted_max_q + 1)) {
    //
    for (int q = extracted_min_q; q <= extracted_min_q; q++) {
        kets.push_back(std::make_shared<const ket_t>(ket_creator.set_quantum_number_q(q).create()));
    }
    return BasisClassicalLight<Scalar>(std::move(kets));
}

// Explicit instantiations
template class BasisClassicalLightCreator<float>;
template class BasisClassicalLightCreator<double>;
template class BasisClassicalLightCreator<std::complex<float>>;
template class BasisClassicalLightCreator<std::complex<double>>;

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include "ket/KetClassicalLightCreator.hpp"
#include "utils/streamed.hpp"
#include <doctest/doctest.h>
#include <spdlog/spdlog.h>

DOCTEST_TEST_CASE("create a classical light basis") {
    float energy{0.03};
    auto basis = BasisClassicalLightCreator<float>()
                     .set_photon_energy(energy)
                     .restrict_quantum_number_q(-3, 3)
                     .create();
    for (const auto &ket : basis) {
        DOCTEST_CHECK(ket.get_photon_energy() == energy);
        DOCTEST_CHECK(ket.get_energy() == ket.get_photon_energy() * ket.get_quantum_number_q());
        SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Ket: {}", fmt::streamed(ket));
    }
}
