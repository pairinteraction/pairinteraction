// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <optional>
#include <type_traits>

namespace pairinteraction {
struct AtomDescriptionByParameters {
    Parity parity{Parity::UNKNOWN};
    std::optional<double> energy;
    std::optional<double> quantum_number_f;
    std::optional<double> quantum_number_m;
    std::optional<int> quantum_number_n;
    std::optional<double> quantum_number_nu;
    std::optional<double> quantum_number_nui;
    std::optional<double> quantum_number_l;
    std::optional<double> quantum_number_s;
    std::optional<double> quantum_number_j;
    std::optional<double> quantum_number_l_ryd;
    std::optional<double> quantum_number_j_ryd;
};
} // namespace pairinteraction
