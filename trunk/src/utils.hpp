#ifndef __UTILS_HPP
#define __UTILS_HPP

namespace utils {

template<typename T>
std::vector<T> linspace(T a, T b, int n) {
  std::vector<T> array(n);
  T step = (b - a) / (n - 1);
  for (int i = 0; i != n; ++i)
    array[i] = a + i*step;
  return array;
}

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
