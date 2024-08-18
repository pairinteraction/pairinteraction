#pragma once

#include "pintr/enums/Parity.hpp"
#include "pintr/utils/Range.hpp"
#include "pintr/utils/traits.hpp"

#include <optional>
#include <type_traits>

namespace pintr {
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
} // namespace pintr
