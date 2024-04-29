#pragma once

enum class SortBy : unsigned char {
    NONE = 1 << 0,
    QUANTUM_NUMBER_F = 1 << 1,
    QUANTUM_NUMBER_M = 1 << 2,
    PARITY = 1 << 3,
    KET = 1 << 4,
    DIAGONAL = 1 << 5
};

inline constexpr SortBy operator&(SortBy x, SortBy y) {
    return static_cast<SortBy>(static_cast<unsigned char>(x) & static_cast<unsigned char>(y));
}

inline constexpr SortBy operator|(SortBy x, SortBy y) {
    return static_cast<SortBy>(static_cast<unsigned char>(x) | static_cast<unsigned char>(y));
}

inline constexpr SortBy operator~(SortBy x) {
    return static_cast<SortBy>(~static_cast<unsigned char>(x));
}
