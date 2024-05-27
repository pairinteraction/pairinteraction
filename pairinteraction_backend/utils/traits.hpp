#pragma once

#include <complex>
#include <type_traits>

namespace traits {

/**
 * @struct CrtpTraits
 *
 * @brief Helper struct to extract types from a derived basis type. Must be specialized for each
 * derived basis type.
 *
 * @tparam Derived The derived basis type from which to extract types.
 */

template <typename Derived>
struct CrtpTraits;

/**
 * @struct NumTraits
 *
 * @brief Helper struct to extract types from a numerical type.
 *
 * @tparam Numeric The numerical type from which to extract types.
 */

template <typename Numeric>
struct NumTraits {
    static_assert(std::is_floating_point_v<Numeric>);
    using real_t = Numeric;
    static constexpr bool is_complex = false;
};

template <typename Numeric>
struct NumTraits<std::complex<Numeric>> {
    static_assert(std::is_floating_point_v<Numeric>);
    using real_t = Numeric;
    static constexpr bool is_complex = true;
};

} // namespace traits
