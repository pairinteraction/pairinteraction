// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <optional>
#include <type_traits>

namespace pairinteraction {
struct AtomDescriptionByRanges {
    Parity parity{Parity::UNKNOWN};
    Range<double> range_energy;
    Range<double> range_quantum_number_f;
    Range<double> range_quantum_number_m;
    Range<int> range_quantum_number_n;
    Range<double> range_quantum_number_nu;
    Range<double> range_quantum_number_nui;
    Range<double> range_quantum_number_l;
    Range<double> range_quantum_number_s;
    Range<double> range_quantum_number_j;
    Range<double> range_quantum_number_l_ryd;
    Range<double> range_quantum_number_j_ryd;
};
} // namespace pairinteraction
