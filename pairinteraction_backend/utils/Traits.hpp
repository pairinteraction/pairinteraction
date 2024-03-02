#pragma once

#include <complex>

namespace internal {

/**
 * @struct BasisTraits
 *
 * @brief Helper struct to extract types from a derived basis type. Must be specialized for each
 * derived basis type.
 *
 * @tparam Derived The derived basis type from which to extract types.
 */

template <typename Derived>
struct BasisTraits;

/**
 * @struct NumTraits
 *
 * @brief Helper struct to extract types from a numerical type.
 *
 * @tparam Numeric The numerical type from which to extract types.
 */

template <typename Numeric>
struct NumTraits {
    using real_t = Numeric;
};

template <typename Numeric>
struct NumTraits<std::complex<Numeric>> {
    using real_t = Numeric;
};

} // namespace internal
