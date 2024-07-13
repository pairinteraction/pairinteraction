#pragma once

#include <functional>
#include <initializer_list>

namespace pintr {
template <typename Derived>
class System;

template <typename Scalar>
class DiagonalizerInterface;

template <typename Derived>
void diagonalize(std::initializer_list<std::reference_wrapper<System<Derived>>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 int precision = 4);

template <typename Derived>
void diagonalize(std::initializer_list<std::reference_wrapper<System<Derived>>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 typename Derived::real_t min_eigenvalue, typename Derived::real_t max_eigenvalue,
                 int precision = 4);
} // namespace pintr
