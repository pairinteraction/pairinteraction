#pragma once

#include <functional>
#include <initializer_list>
#include <vector>

namespace pintr {
template <typename Scalar>
class DiagonalizerInterface;

template <typename Derived>
void diagonalize(std::initializer_list<std::reference_wrapper<Derived>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 int precision = 12);

template <typename Derived>
void diagonalize(std::initializer_list<std::reference_wrapper<Derived>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 typename Derived::real_t min_eigenvalue, typename Derived::real_t max_eigenvalue,
                 int precision = 12);

template <typename Derived>
void diagonalize(std::vector<Derived> &systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 int precision = 12);

template <typename Derived>
void diagonalize(std::vector<Derived> &systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 typename Derived::real_t min_eigenvalue, typename Derived::real_t max_eigenvalue,
                 int precision = 12);
} // namespace pintr
