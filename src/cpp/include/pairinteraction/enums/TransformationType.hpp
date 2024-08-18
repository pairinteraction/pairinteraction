#pragma once

#include <initializer_list>
#include <vector>

namespace pintr {
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
inline TransformationType &operator|=(TransformationType &x, TransformationType y) {
    return x = x | y;
}
inline TransformationType &operator&=(TransformationType &x, TransformationType y) {
    return x = x & y;
}

namespace utils {
inline bool is_sorting(TransformationType label) {
    static constexpr TransformationType MASK_SORTING = TransformationType::SORT_BY_KET |
        TransformationType::SORT_BY_QUANTUM_NUMBER_F |
        TransformationType::SORT_BY_QUANTUM_NUMBER_M | TransformationType::SORT_BY_PARITY |
        TransformationType::SORT_BY_ENERGY;
    return (label & ~MASK_SORTING) == TransformationType::IDENTITY;
}

inline bool has_bit(TransformationType value, TransformationType bit) {
    return (value & bit) == bit;
}

inline bool is_comprised_by_label(TransformationType label,
                                  const std::vector<TransformationType> &list_used) {
    TransformationType used = TransformationType::IDENTITY;
    for (auto l : list_used) {
        used |= l;
    }
    return label == used;
}

inline bool is_sorted_by_label(TransformationType label,
                               const std::vector<TransformationType> &list_used) {
    TransformationType used = TransformationType::IDENTITY;
    for (auto it = list_used.rbegin(); it != list_used.rend(); ++it) {
        if ((label & *it) != *it) {
            break;
        }
        used |= *it;
    }
    return label == used;
}

} // namespace utils
} // namespace pintr
