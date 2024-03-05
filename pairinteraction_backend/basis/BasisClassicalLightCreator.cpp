#include "basis/BasisClassicalLightCreator.hpp"

template <typename Scalar>
BasisClassicalLightCreator<Scalar>::BasisClassicalLightCreator() {
} // TODO: formal empty constructor

template <typename Scalar>
BasisClassicalLightCreator<Scalar> &
BasisClassicalLightCreator<Scalar>::restrict_quantum_number_n_sigma_p(int min, int max) {
    min_quantum_number_n_sigma_p.emplace(min);
    max_quantum_number_n_sigma_p.emplace(max);
    return *this;
}

template <typename Scalar>
BasisClassicalLightCreator<Scalar> &
BasisClassicalLightCreator<Scalar>::restrict_quantum_number_n_pi(int min, int max) {
    min_quantum_number_n_pi.emplace(min);
    max_quantum_number_n_pi.emplace(max);
    return *this;
}

template <typename Scalar>
BasisClassicalLightCreator<Scalar> &
BasisClassicalLightCreator<Scalar>::restrict_quantum_number_n_sigma_m(int min, int max) {
    min_quantum_number_n_sigma_m.emplace(min);
    max_quantum_number_n_sigma_m.emplace(max);
    return *this;
}

template <typename Scalar>
BasisClassicalLight<Scalar> BasisClassicalLightCreator<Scalar>::create() const {
    // TODO: Proper implementation

    std::vector<std::shared_ptr<const ket_t>> kets;
    kets.reserve(2);
    kets.push_back(
        std::make_shared<const ket_t>(KetClassicalLightCreator<real_t>(20, 1, 0, 0).create()));
    kets.push_back(
        std::make_shared<const ket_t>(KetClassicalLightCreator<real_t>(20, 0, 0, 0).create()));
    return BasisClassicalLight<Scalar>(std::move(kets));
}

// Explicit instantiations
template class BasisClassicalLightCreator<float>;
template class BasisClassicalLightCreator<double>;
template class BasisClassicalLightCreator<std::complex<float>>;
template class BasisClassicalLightCreator<std::complex<double>>;
