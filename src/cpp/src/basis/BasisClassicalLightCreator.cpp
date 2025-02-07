#include "pairinteraction/basis/BasisClassicalLightCreator.hpp"

#include "pairinteraction/basis/BasisClassicalLight.hpp"
#include "pairinteraction/ket/KetClassicalLightCreator.hpp"

#include <stdexcept>

namespace pairinteraction {
template <typename Scalar>
BasisClassicalLightCreator<Scalar> &
BasisClassicalLightCreator<Scalar>::set_photon_energy(real_t value) {

    if (value < 0) {
        throw std::invalid_argument("photon energy must be >= 0");
    }

    photon_energy = value;
    return *this;
}

template <typename Scalar>
BasisClassicalLightCreator<Scalar> &
BasisClassicalLightCreator<Scalar>::restrict_quantum_number_q(int min, int max) {

    if (min > max) {
        throw std::invalid_argument("min has to be smaller then max");
    }
    range_quantum_number_q = {min, max};
    return *this;
}

template <typename Scalar>
std::shared_ptr<const BasisClassicalLight<Scalar>>
BasisClassicalLightCreator<Scalar>::create() const {

    std::vector<std::shared_ptr<const ket_t>> kets;
    kets.reserve(range_quantum_number_q.max() - range_quantum_number_q.min() + 1);

    auto ket_creator = KetClassicalLightCreator<real_t>().set_photon_energy(photon_energy);

    for (int q = range_quantum_number_q.min(); q <= range_quantum_number_q.max(); q++) {
        kets.push_back(ket_creator.set_quantum_number_q(q).create());
    }

    return std::make_shared<const BasisClassicalLight<Scalar>>(
        typename BasisClassicalLight<Scalar>::Private(), std::move(kets));
}

// Explicit instantiations
template class BasisClassicalLightCreator<double>;
template class BasisClassicalLightCreator<std::complex<double>>;
} // namespace pairinteraction
