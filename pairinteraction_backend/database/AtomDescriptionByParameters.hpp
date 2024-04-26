#pragma once

#include "utils/traits.hpp"

#include <optional>

template <typename Scalar>
struct AtomDescriptionByParameters {
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    std::optional<real_t> energy;
    std::optional<float> quantum_number_f;
    std::optional<float> quantum_number_m;
    std::optional<int> parity;
    std::optional<int> quantum_number_n;
    std::optional<real_t> quantum_number_nu;
    std::optional<real_t> quantum_number_l;
    std::optional<real_t> quantum_number_s;
    std::optional<real_t> quantum_number_j;
};
