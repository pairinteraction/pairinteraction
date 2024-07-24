#pragma once

#include <complex>
#include <functional>
#include <type_traits>

namespace pintr::traits {

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

/**
 * @struct OpTraits
 *
 * @brief Helper struct to extract whether a type supports certain operations.
 *
 * @tparam T The type for which the existence of operations is to be checked.
 */

template <typename T>
struct OpTraits {
    static constexpr bool has_equal_v = std::is_invocable_r_v<bool, std::equal_to<>, T, T>;
    static constexpr bool has_less_v = std::is_invocable_r_v<bool, std::less<>, T, T>;
};

} // namespace pintr::traits
