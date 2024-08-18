#pragma once

#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <optional>
#include <type_traits>

namespace pairinteraction {
template <typename Real>
struct AtomDescriptionByParameters {
    static_assert(std::is_floating_point_v<Real>);

    Parity parity{Parity::UNKNOWN};
    std::optional<Real> energy;
    std::optional<Real> quantum_number_f;
    std::optional<Real> quantum_number_m;
    std::optional<int> quantum_number_n;
    std::optional<Real> quantum_number_nu;
    std::optional<Real> quantum_number_l;
    std::optional<Real> quantum_number_s;
    std::optional<Real> quantum_number_j;
};
} // namespace pairinteraction
