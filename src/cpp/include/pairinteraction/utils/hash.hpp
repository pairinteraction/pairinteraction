// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <array>
#include <complex>
#include <functional>
#include <vector>

namespace pairinteraction::utils {

/**
 * @struct hash
 *
 * @brief Hash function
 *
 * The `std::hash` template allows specialization but only for types that are
 * not in the standard library. This means that we cannot specialize
 * `std::hash` for, e.g. `std::array`. To this end we define a struct `hash`
 * which just inherits from `std::hash` by default.
 *
 * @tparam T type to be hashed
 */

template <typename T>
struct hash;

/**
 * @function hash_combine
 *
 * @brief Combine hashes
 *
 * The implementation of `hash_combine` is copied from Boost but simplified.
 *
 * @param seed  start hash
 * @param v  value whose hash is to be added to \p seed
 *
 * @tparam T type to be hashed
 */

template <typename T>
inline void hash_combine(std::size_t &seed, T const &v) {
    hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/**
 * @function hash_range
 *
 * @brief Combine hashes of values in a range
 *
 * @param first forward iterator
 * @param last forward iterator
 *
 * @return combined hash
 *
 * @tparam It forward iterator
 */

template <typename It>
inline std::size_t hash_range(It first, It last) {
    std::size_t seed = 0;
    for (; first != last; ++first) {
        hash_combine(seed, *first);
    }
    return seed;
}

// By default use std::hash
template <typename T>
struct hash : std::hash<T> {};

// Specializations for other types
template <typename T, std::size_t N>
struct hash<std::array<T, N>> {
    std::size_t operator()(std::array<T, N> const &a) const {
        return hash_range(a.begin(), a.end());
    }
};

template <typename T>
struct hash<std::vector<T>> {
    std::size_t operator()(std::vector<T> const &v) const { return hash_range(v.begin(), v.end()); }
};

template <typename T>
struct hash<std::complex<T>> {
    std::size_t operator()(std::complex<T> const &c) const {
        std::size_t seed = 0;
        hash_combine(seed, c.real());
        hash_combine(seed, c.imag());
        return seed;
    }
};

enum class Parity : int;

template <>
struct hash<Parity> {
    std::size_t operator()(const Parity &parity) const {
        return std::hash<char>{}(static_cast<char>(parity));
    }
};

} // namespace pairinteraction::utils
