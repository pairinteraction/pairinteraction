// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/enums/FloatType.hpp"

#include <complex>
#include <functional>
#include <type_traits>

namespace pairinteraction::traits {

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

/**
 * @struct FPTypeTraits
 *
 * @brief Helper struct to extract the type from a floating point precision enum.
 *
 * @tparam T The enum value for which the type is to be extracted.
 */

template <FloatType T>
struct FPTypeTraits;

template <>
struct FPTypeTraits<FloatType::FLOAT32> {
    using type = float;
};

template <>
struct FPTypeTraits<FloatType::FLOAT64> {
    using type = double;
};

template <typename Scalar, FloatType FP>
using restricted_t = std::conditional_t<traits::NumTraits<Scalar>::is_complex_v,
                                        std::complex<typename FPTypeTraits<FP>::type>,
                                        typename FPTypeTraits<FP>::type>;

} // namespace pairinteraction::traits
