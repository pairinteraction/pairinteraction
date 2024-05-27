#pragma once

#include <initializer_list>
#include <vector>

enum class TransformationType : unsigned char {
    NONE = 0,
    SORT_BY_KET = 1 << 0,
    SORT_BY_QUANTUM_NUMBER_F = 1 << 1,
    SORT_BY_QUANTUM_NUMBER_M = 1 << 2,
    SORT_BY_PARITY = 1 << 3,
    SORT_BY_ENERGY = 1 << 4,
    ROTATE = 1 << 5,
    ARBITRARY = 1 << 6,
    MASK_SORTING = SORT_BY_KET | SORT_BY_QUANTUM_NUMBER_F | SORT_BY_QUANTUM_NUMBER_M |
        SORT_BY_PARITY | SORT_BY_ENERGY
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
    return (label & ~TransformationType::MASK_SORTING) == TransformationType::NONE;
}

inline bool has_bit(TransformationType value, TransformationType bit) {
    return (value & bit) == bit;
}

inline bool is_comprised_by_label(TransformationType label,
                                  const std::vector<TransformationType> &list_used) {
    TransformationType used = TransformationType::NONE;
    for (auto l : list_used) {
        used |= l;
    }
    return label == used;
}

inline bool is_sorted_by_label(TransformationType label,
                               const std::vector<TransformationType> &list_used) {
    TransformationType used = TransformationType::NONE;
    for (auto it = list_used.rbegin(); it != list_used.rend(); ++it) {
        if ((label & *it) != *it) {
            break;
        }
        used |= *it;
    }
    return label == used;
}

} // namespace utils
