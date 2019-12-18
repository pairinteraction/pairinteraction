/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * This file is part of the pairinteraction library.
 *
 * The pairinteraction library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The pairinteraction library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <array>
#include <complex>
#include <functional>
#include <random>
#include <type_traits>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#else
#include <unistd.h>
#endif

namespace utils {

template <typename T>
struct is_complex : std::false_type {};
template <typename S>
struct is_complex<std::complex<S>> : std::true_type {};

template <typename T>
inline T conjugate(const T &val) {
    return val;
}
template <typename S>
inline std::complex<S> conjugate(const std::complex<S> &val) {
    return std::conj(val);
}

template <typename T>
inline bool is_true(const T &val) {
    return std::all_of(val.begin(), val.end(), [](bool i) { return i; });
}
inline bool is_true(bool val) { return val; }

template <typename T>
inline typename std::enable_if<!is_complex<T>::value, T>::type imaginary_unit() {
    throw std::runtime_error(
        "For operations that invoke the imaginary number, a complex data type is needed.");
}
template <typename T>
inline typename std::enable_if<is_complex<T>::value, T>::type imaginary_unit() {
    return {0, 1};
}

template <typename T, typename U>
inline T convert(const U &val) {
    return val;
}
template <typename S>
inline S convert(const std::complex<S> &val) {
    return std::real(val);
}

/** \brief Thread-local static random engine
 *
 * To save some effort this function initializes a static thread-local
 * random engine to avoid race conditions and broken random sampling.
 * It is seeded once using `std::random_device` for good entropy.
 *
 * \returns Reference to static thread-local random engine
 */
inline std::default_random_engine &randint_engine() {
    static thread_local std::default_random_engine eng{std::random_device{}()};
    return eng;
}

/** \brief Generate a random integer
 *
 * This is very similar to the implementation of randint in the GCC
 * standard library Fundamentals TS v2.  It is specified by the
 * standard per clause 13.2.2.1, Function template randint.
 *
 * The function generates a random integer in the closed interval
 * [\p a,\p b].
 *
 * \param a  lower bound
 * \param b  upper bound
 * \returns random integer between \p a and \p b.
 */
template <typename T>
inline T randint(T a, T b) {
    static_assert(std::is_integral<T>::value && sizeof(T) > 1, "The type must be an integer!");
    return std::uniform_int_distribution<T>(a, b)(randint_engine());
}

// https://de.wikipedia.org/wiki/FNV_(Informatik)
inline uint64_t FNV64(const uint8_t *s, size_t sz) {
    const uint64_t magicPrime = 0x00000100000001b3;
    uint64_t hash = 0xcbf29ce484222325;

    for (size_t i = 0; i < sz; ++i) {
        hash = (hash ^ s[i]) * magicPrime;
    }
    return hash;
}

inline long get_pid() {
#ifdef _WIN32
    return GetCurrentProcessId();
#else
    return ::getpid();
#endif
}

/// \brief Hash function
///
/// The `std::hash` template allows specialization but only for types that are
/// not in the standard library.  This means that we cannot specialize
/// `std::hash` for, e.g. `std::array`.  To this end we define a struct `hash`
/// which just inherits from `std::hash` by default.
///
/// \tparam Key  type to be hashed
template <typename Key>
struct hash;

/// \brief Combine hashes
///
/// The implementation of `hash_combine` is copied from Boost but simplified.
/// It uses the custom `hash`.
///
/// \param seed  start hash
/// \param v  value whose hash is to be added to \p seed
template <typename T>
inline void hash_combine(std::size_t &seed, T const &v) {
    hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/// \brief Combine hashes of values in a range
///
/// \param first  forward iterator
/// \param last  forward iterator
/// \returns combined hash of all values in the range
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
    using argument_type = std::array<T, N>;
    using result_type = std::size_t;
    std::size_t operator()(std::array<T, N> const &a) const {
        return hash_range(a.begin(), a.end());
    }
};

template <typename T>
struct hash<std::complex<T>> {
    using argument_type = std::complex<T>;
    using result_type = std::size_t;
    std::size_t operator()(std::complex<T> const &c) const {
        std::size_t seed = 0;
        hash_combine(seed, c.real());
        hash_combine(seed, c.imag());
        return seed;
    }
};

} // namespace utils

#endif // UTILS_H
