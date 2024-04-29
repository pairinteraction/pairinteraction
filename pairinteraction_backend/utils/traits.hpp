#pragma once

#include <complex>

namespace traits {

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
 * @struct OperatorTraits
 *
 * @brief Helper struct to extract types from a derived operator type. Must be specialized for each
 * derived operator type.
 *
 * @tparam Derived The derived operator type from which to extract types.
 */

template <typename Derived>
struct OperatorTraits;

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

/**
 * @struct is_complex
 *
 * @brief Helper struct to check whether a scalar type is complex.
 *
 * @tparam Scalar The scalar type to check.
 */

template <typename Scalar>
struct is_complex : public std::false_type {};

template <typename Scalar>
struct is_complex<std::complex<Scalar>> : public std::true_type {};

template <typename Scalar>
inline constexpr bool is_complex_v = is_complex<Scalar>::value;

} // namespace traits
