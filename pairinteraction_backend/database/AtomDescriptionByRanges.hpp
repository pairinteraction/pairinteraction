#pragma once

#include "utils/traits.hpp"

#include <optional>
#include <type_traits>

template <typename Real>
struct AtomDescriptionByRanges {
    static_assert(std::is_floating_point_v<Real>);

    std::optional<Real> min_energy;
    std::optional<Real> max_energy;
    std::optional<Real> min_quantum_number_f;
    std::optional<Real> max_quantum_number_f;
    std::optional<Real> min_quantum_number_m;
    std::optional<Real> max_quantum_number_m;
    std::optional<int> parity;
    std::optional<int> min_quantum_number_n;
    std::optional<int> max_quantum_number_n;
    std::optional<Real> min_quantum_number_nu;
    std::optional<Real> max_quantum_number_nu;
    std::optional<Real> min_quantum_number_l;
    std::optional<Real> max_quantum_number_l;
    std::optional<Real> min_quantum_number_s;
    std::optional<Real> max_quantum_number_s;
    std::optional<Real> min_quantum_number_j;
    std::optional<Real> max_quantum_number_j;
};
