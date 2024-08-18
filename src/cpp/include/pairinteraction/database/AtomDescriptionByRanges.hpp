#pragma once

#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <optional>
#include <type_traits>

namespace pairinteraction {
template <typename Real>
struct AtomDescriptionByRanges {
    static_assert(std::is_floating_point_v<Real>);

    Parity parity{Parity::UNKNOWN};
    Range<Real> range_energy;
    Range<Real> range_quantum_number_f;
    Range<Real> range_quantum_number_m;
    Range<int> range_quantum_number_n;
    Range<Real> range_quantum_number_nu;
    Range<Real> range_quantum_number_l;
    Range<Real> range_quantum_number_s;
    Range<Real> range_quantum_number_j;
};
} // namespace pairinteraction
