#pragma once

#include "utils/traits.hpp"

#include <optional>

template <typename Scalar>
struct AtomDescriptionByRanges {
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    std::optional<real_t> min_energy;
    std::optional<real_t> max_energy;
    std::optional<float> min_quantum_number_f;
    std::optional<float> max_quantum_number_f;
    std::optional<float> min_quantum_number_m;
    std::optional<float> max_quantum_number_m;
    std::optional<int> parity;
    std::optional<int> min_quantum_number_n;
    std::optional<int> max_quantum_number_n;
    std::optional<real_t> min_quantum_number_nu;
    std::optional<real_t> max_quantum_number_nu;
    std::optional<real_t> min_quantum_number_l;
    std::optional<real_t> max_quantum_number_l;
    std::optional<real_t> min_quantum_number_s;
    std::optional<real_t> max_quantum_number_s;
    std::optional<real_t> min_quantum_number_j;
    std::optional<real_t> max_quantum_number_j;
};
