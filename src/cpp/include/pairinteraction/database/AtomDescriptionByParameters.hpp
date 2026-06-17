// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/traits.hpp"

#include <optional>
#include <string>
#include <type_traits>
#include <unordered_map>

namespace pairinteraction {
struct AtomDescriptionByParameters {
    std::optional<double> energy;
    // Quantum numbers keyed by their logical name (e.g. "f", "m", "n", "nu", "l", ..., "parity").
    std::unordered_map<std::string, double> quantum_numbers;
};
} // namespace pairinteraction
