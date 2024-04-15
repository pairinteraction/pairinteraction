#pragma once

#include "basis/BasisClassicalLight.hpp"
#include "ket/KetClassicalLight.hpp"
#include "utils/Traits.hpp"

#include <complex>
#include <limits>
#include <optional>
#include <string>

template <typename Scalar>
class BasisClassicalLightCreator {
public:
    using real_t = typename Traits::NumTraits<Scalar>::real_t;
    using ket_t = KetClassicalLight<real_t>;
    BasisClassicalLightCreator();
    BasisClassicalLightCreator<Scalar> &set_energy(Scalar energy);
    BasisClassicalLightCreator<Scalar> &restrict_quantum_number_n_sigma_p(int min, int max);
    BasisClassicalLightCreator<Scalar> &restrict_quantum_number_n_pi(int min, int max);
    BasisClassicalLightCreator<Scalar> &restrict_quantum_number_n_sigma_m(int min, int max);
    BasisClassicalLight<Scalar> create() const;

private:
    std::optional<Scalar> energy;
    std::optional<int> min_quantum_number_n_sigma_p;
    std::optional<int> max_quantum_number_n_sigma_p;
    std::optional<int> min_quantum_number_n_pi;
    std::optional<int> max_quantum_number_n_pi;
    std::optional<int> min_quantum_number_n_sigma_m;
    std::optional<int> max_quantum_number_n_sigma_m;
};

extern template class BasisClassicalLightCreator<float>;
extern template class BasisClassicalLightCreator<double>;
extern template class BasisClassicalLightCreator<std::complex<float>>;
extern template class BasisClassicalLightCreator<std::complex<double>>;
