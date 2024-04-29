#include "basis/BasisClassicalLightCreator.hpp"

#include "basis/BasisClassicalLight.hpp"
#include "ket/KetClassicalLightCreator.hpp"

template <typename Scalar>
BasisClassicalLightCreator<Scalar> &
BasisClassicalLightCreator<Scalar>::restrict_quantum_number_q(int min, int max) {
    min_quantum_number_q.emplace(min);
    max_quantum_number_q.emplace(max);
    return *this;
}

template <typename Scalar>
BasisClassicalLight<Scalar> BasisClassicalLightCreator<Scalar>::create() const {
    // TODO: Proper implementation

    std::vector<std::shared_ptr<const ket_t>> kets;
    kets.reserve(2);
    kets.push_back(std::make_shared<const ket_t>(KetClassicalLightCreator<real_t>(20, 1).create()));
    kets.push_back(std::make_shared<const ket_t>(KetClassicalLightCreator<real_t>(20, 0).create()));
    return BasisClassicalLight<Scalar>(std::move(kets));
}

// Explicit instantiations
template class BasisClassicalLightCreator<float>;
template class BasisClassicalLightCreator<double>;
template class BasisClassicalLightCreator<std::complex<float>>;
template class BasisClassicalLightCreator<std::complex<double>>;
