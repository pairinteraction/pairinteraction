#pragma once

#include "utils/Range.hpp"
#include "utils/traits.hpp"

#include <complex>
#include <limits>
#include <memory>
#include <optional>
#include <ranges>
#include <string>
#include <vector>

template <typename Scalar>
class BasisClassicalLight;

template <typename Real>
class KetClassicalLight;

template <typename Scalar>
class BasisClassicalLightCreator {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

public:
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using ket_t = KetClassicalLight<real_t>;
    BasisClassicalLightCreator() = default;
    BasisClassicalLightCreator<Scalar> &set_photon_energy(real_t value);
    BasisClassicalLightCreator<Scalar> &restrict_quantum_number_q(int min, int max);
    std::shared_ptr<const BasisClassicalLight<Scalar>> create() const;

private:
    friend class BasisClassicalLight<Scalar>;
    // friend class KetClassicalLight<Scalar>;
    real_t photon_energy{0.0};
    Range<int> range_quantum_number_q{0, 0};
};

extern template class BasisClassicalLightCreator<float>;
extern template class BasisClassicalLightCreator<double>;
extern template class BasisClassicalLightCreator<std::complex<float>>;
extern template class BasisClassicalLightCreator<std::complex<double>>;
