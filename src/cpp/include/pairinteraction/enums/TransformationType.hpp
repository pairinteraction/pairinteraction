// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <set>
#include <vector>

namespace pairinteraction {
enum class TransformationType : unsigned char {
    IDENTITY = 0,
    SORT_BY_KET = 1 << 0,
    SORT_BY_QUANTUM_NUMBER_F = 1 << 1,
    SORT_BY_QUANTUM_NUMBER_M = 1 << 2,
    SORT_BY_PARITY = 1 << 3,
    SORT_BY_ENERGY = 1 << 4,
    ROTATE = 1 << 5,
    ARBITRARY = 1 << 6
};

inline constexpr TransformationType operator&(TransformationType x, TransformationType y) {
    return static_cast<TransformationType>(static_cast<unsigned char>(x) &
                                           static_cast<unsigned char>(y));
}
inline constexpr TransformationType operator|(TransformationType x, TransformationType y) {
    return static_cast<TransformationType>(static_cast<unsigned char>(x) |
                                           static_cast<unsigned char>(y));
}
inline constexpr TransformationType operator~(TransformationType x) {
    return static_cast<TransformationType>(~static_cast<unsigned char>(x));
}

namespace utils {
inline bool is_sorting(TransformationType label) {
    static constexpr TransformationType MASK_SORTING = TransformationType::SORT_BY_KET |
        TransformationType::SORT_BY_QUANTUM_NUMBER_F |
        TransformationType::SORT_BY_QUANTUM_NUMBER_M | TransformationType::SORT_BY_PARITY |
        TransformationType::SORT_BY_ENERGY;
    return (label & ~MASK_SORTING) == TransformationType::IDENTITY;
}

} // namespace utils
} // namespace pairinteraction
