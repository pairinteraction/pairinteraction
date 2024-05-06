#pragma once

#include <complex>
#include <limits>
#include <optional>
#include <ranges>
#include <string>
#include <vector>

#include "utils/traits.hpp"

template <typename Scalar>
class BasisClassicalLight;

template <typename Real>
class KetClassicalLight;

template <typename Scalar>
class BasisClassicalLightCreator {
public:
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetClassicalLight<real_t>;
    BasisClassicalLightCreator() = default;
    BasisClassicalLightCreator<Scalar> &set_photon_energy(real_t photon_energy);
    BasisClassicalLightCreator<Scalar> &restrict_quantum_number_q(int min, int max);
    BasisClassicalLightCreator<Scalar> &add_ket(const ket_t &ket);
    BasisClassicalLight<Scalar> create() const;

private:
    friend class BasisClassicalLight<Scalar>;
    // friend class KetClassicalLight<Scalar>;
    std::optional<real_t> photon_energy;
    std::optional<int> min_quantum_number_q;
    std::optional<int> max_quantum_number_q;
};

extern template class BasisClassicalLightCreator<float>;
extern template class BasisClassicalLightCreator<double>;
extern template class BasisClassicalLightCreator<std::complex<float>>;
extern template class BasisClassicalLightCreator<std::complex<double>>;
