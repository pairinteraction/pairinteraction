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

}

#endif
