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
    static constexpr bool is_complex_v = false;
};

template <typename Numeric>
struct NumTraits<std::complex<Numeric>> {
    static_assert(std::is_floating_point_v<Numeric>);

    using real_t = Numeric;
    static constexpr bool is_complex_v = true;
};

/**
 * @struct is_complex_or_floating_point
 *
 * @brief Helper struct to check whether a type is a complex or floating point type.
 *
 * @tparam Scalar The scalar type to check.
 */

template <class Scalar>
struct is_complex_or_floating_point : std::is_floating_point<Scalar> {};

template <class Scalar>
struct is_complex_or_floating_point<std::complex<Scalar>> : std::is_floating_point<Scalar> {};

template <typename Scalar>
inline constexpr bool is_complex_or_floating_point_v = is_complex_or_floating_point<Scalar>::value;

} // namespace traits
