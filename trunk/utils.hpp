#ifndef __UTILS_HPP
#define __UTILS_HPP

#include <complex>
#include <type_traits>
#include <functional>

namespace utils {

template<typename T> struct is_complex : public std::false_type {};
template<typename T> struct is_complex<std::complex<T>> : public std::true_type {};

// https://de.wikipedia.org/wiki/FNV_(Informatik)
inline uint64_t FNV64(uint8_t *s, size_t sz)
{
    const uint64_t magicPrime = 0x00000100000001b3;
    uint64_t hash = 0xcbf29ce484222325;

    for(size_t i = 0; i < sz; ++i)
    {
        hash = (hash ^ s[i]) * magicPrime;
    }
    return hash;
}

// from boost
template <class T>
inline void hash_combine(size_t & seed, const T & v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

}

#endif
