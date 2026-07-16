// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <set>
#include <vector>

namespace pairinteraction {
enum class TransformationType : unsigned char {
    CANONICAL_ORDER,
    SORT_BY_QUANTUM_NUMBER_F,
    SORT_BY_QUANTUM_NUMBER_M,
    SORT_BY_PARITY,
    SORT_BY_ENERGY,
    ARBITRARY
};

namespace utils {
inline bool is_sorting(TransformationType label) {
    switch (label) {
    case TransformationType::SORT_BY_QUANTUM_NUMBER_F:
    case TransformationType::SORT_BY_QUANTUM_NUMBER_M:
    case TransformationType::SORT_BY_PARITY:
    case TransformationType::SORT_BY_ENERGY:
        return true;
    default:
        return false;
    }
}

} // namespace utils
} // namespace pairinteraction
