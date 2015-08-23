#ifndef __UTILS_HPP
#define __UTILS_HPP

namespace utils {

// https://de.wikipedia.org/wiki/FNV_(Informatik)
uint64_t FNV64(byte_t *s, size_t sz)
{
    const uint64_t magicPrime = 0x00000100000001b3;
    uint64_t hash = 0xcbf29ce484222325;

    for(size_t i = 0; i < sz; ++i)
    {
        hash = (hash ^ s[i]) * magicPrime;
    }
    return hash;
}

}

#endif
