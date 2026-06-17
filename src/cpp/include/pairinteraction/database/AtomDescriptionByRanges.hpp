// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/Range.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <optional>
#include <string>
#include <type_traits>
#include <unordered_map>

namespace pairinteraction {
struct AtomDescriptionByRanges {
    Range<double> range_energy;
    // Quantum number ranges keyed by their logical name (e.g. "f", "m", "n", "l", ..., "parity").
    std::unordered_map<std::string, Range<double>> quantum_number_ranges;
    double quantum_number_standard_deviation_factor{2};
};
} // namespace pairinteraction
