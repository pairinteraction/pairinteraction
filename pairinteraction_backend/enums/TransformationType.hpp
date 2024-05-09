#pragma once

#include <initializer_list>

enum class TransformationType : unsigned char {
    NONE = 0,
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
inline TransformationType &operator|=(TransformationType &x, TransformationType y) {
    return x = x | y;
}
inline TransformationType &operator&=(TransformationType &x, TransformationType y) {
    return x = x & y;
}

namespace utils {
inline bool is_sorting(TransformationType label) {
    TransformationType sorting_label = TransformationType::NONE;
    for (auto l : {TransformationType::SORT_BY_QUANTUM_NUMBER_F,
                   TransformationType::SORT_BY_QUANTUM_NUMBER_M, TransformationType::SORT_BY_PARITY,
                   TransformationType::SORT_BY_KET, TransformationType::SORT_BY_ENERGY}) {
        if ((label & l) == l) {
            sorting_label |= l;
        }
    }
    return sorting_label == label;
}
} // namespace utils
