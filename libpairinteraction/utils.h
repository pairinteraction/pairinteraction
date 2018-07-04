/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef UTILS_H
#define UTILS_H

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

namespace utils
{

template <typename T>
struct is_complex : public std::false_type {
};
template <typename T>
struct is_complex<std::complex<T>> : public std::true_type {
};

/** \brief Thread-local static random engine
 *
 * To save some effort this function initializes a static thread-local
 * random engine to avoid race conditions and broken random sampling.
 * It is seeded once using `std::random_device` for good entropy.
 *
 * \returns Reference to static thread-local random engine
 */
inline std::default_random_engine &randint_engine()
{
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
inline T randint(T a, T b)
{
    static_assert(std::is_integral<T>::value && sizeof(T) > 1,
                  "The type must be an integer!");
    return std::uniform_int_distribution<T>(a, b)(randint_engine());
}

// https://de.wikipedia.org/wiki/FNV_(Informatik)
inline uint64_t FNV64(const uint8_t *s, size_t sz)
{
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
}

#endif // UTILS_H
