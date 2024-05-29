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
    static_assert(std::is_arithmetic_v<Numeric>);

    using real_t = Numeric;
    static constexpr bool is_complex_v = false;
    static constexpr bool from_floating_point_v = std::is_floating_point<Numeric>::value;
    static constexpr bool from_integral_v = std::is_integral<Numeric>::value;
};

template <typename Numeric>
struct NumTraits<std::complex<Numeric>> {
    static_assert(std::is_arithmetic_v<Numeric>);

    using real_t = Numeric;
    static constexpr bool is_complex_v = true;
    static constexpr bool from_floating_point_v = std::is_floating_point<Numeric>::value;
    static constexpr bool from_integral_v = std::is_integral<Numeric>::value;
};

} // namespace traits
