#pragma once

#include <functional>
#include <initializer_list>
#include <optional>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class DiagonalizerInterface;

// Note that although a vector is also constructible from a list, the overload resolution
// will prefer the initializer list overload because of less conversions required.

template <typename Derived>
void diagonalize(std::initializer_list<std::reference_wrapper<Derived>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 std::optional<typename Derived::real_t> min_eigenvalue = {},
                 std::optional<typename Derived::real_t> max_eigenvalue = {}, int precision = 12);

template <typename Derived>
void diagonalize(std::vector<Derived> &systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 std::optional<typename Derived::real_t> min_eigenvalue = {},
                 std::optional<typename Derived::real_t> max_eigenvalue = {}, int precision = 12);

template <typename Derived>
void diagonalize(std::vector<std::reference_wrapper<Derived>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 std::optional<typename Derived::real_t> min_eigenvalue = {},
                 std::optional<typename Derived::real_t> max_eigenvalue = {}, int precision = 12);

} // namespace pairinteraction
