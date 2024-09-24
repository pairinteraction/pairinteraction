#pragma once

#include <functional>
#include <initializer_list>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class DiagonalizerInterface;

template <typename Sortable>
class Range;

// Note that although a vector is also constructible from a list, the overload resolution
// will prefer the initializer list overload because of less conversions required.

template <typename Derived>
void diagonalize(std::initializer_list<std::reference_wrapper<Derived>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 int precision = 12, const Range<typename Derived::real_t> &eigenvalue_range = {});

template <typename Derived>
void diagonalize(std::vector<Derived> &systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 int precision = 12, const Range<typename Derived::real_t> &eigenvalue_range = {});

template <typename Derived>
void diagonalize(std::vector<std::reference_wrapper<Derived>> systems,
                 const DiagonalizerInterface<typename Derived::scalar_t> &diagonalizer,
                 int precision = 12, const Range<typename Derived::real_t> &eigenvalue_range = {});

} // namespace pairinteraction
