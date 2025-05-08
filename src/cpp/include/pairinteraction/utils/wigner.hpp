// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/maths.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace pairinteraction::wigner {

namespace {

template <typename Real>
inline constexpr Real PI = 3.141592653589793238462643383279502884;

template <typename Real>
inline Real wigner_uppercase_d_matrix_pi_half(Real f, Real m_initial, Real m_final) {
    static_assert(std::is_floating_point_v<Real>);
    assert(2 * f == std::floor(2 * f) && f >= 0);
    assert(2 * m_initial == std::floor(2 * m_initial) && m_initial >= -f && m_initial <= f);
    assert(2 * m_final == std::floor(2 * m_final) && m_final >= -f && m_final <= f);

    Real result = 0;
    for (int k = std::max(0, static_cast<int>(m_final - m_initial));
         k <= f + std::min(m_final, -m_initial); ++k) {
        result += std::pow(-1, k) * maths::binomial_coefficient(f + m_final, static_cast<Real>(k)) *
            maths::binomial_coefficient(f - m_final, k + m_initial - m_final);
    }
    result *= std::pow(-1, m_initial - m_final) / std::pow(2, f) *
        std::sqrt(maths::factorial(f + m_initial) * maths::factorial(f - m_initial) /
                  (maths::factorial(f + m_final) * maths::factorial(f - m_final)));
    return result;
}

} // namespace

template <typename Scalar>
inline Scalar wigner_uppercase_d_matrix(typename traits::NumTraits<Scalar>::real_t f,
                                        typename traits::NumTraits<Scalar>::real_t m_initial,
                                        typename traits::NumTraits<Scalar>::real_t m_final,
                                        typename traits::NumTraits<Scalar>::real_t alpha,
                                        typename traits::NumTraits<Scalar>::real_t beta,
                                        typename traits::NumTraits<Scalar>::real_t gamma) {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);
    assert(2 * f == std::floor(2 * f) && f >= 0);
    assert(2 * m_initial == std::floor(2 * m_initial) && m_initial >= -f && m_initial <= f);
    assert(2 * m_final == std::floor(2 * m_final) && m_final >= -f && m_final <= f);

    auto scale = std::max(std::abs(m_initial), std::abs(m_final)) *
        std::max(std::abs(alpha), std::abs(gamma));
    auto numerical_precision =
        100 * scale * std::numeric_limits<typename traits::NumTraits<Scalar>::real_t>::epsilon();

    if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
        return Scalar(std::cos(-m_initial * alpha), std::sin(-m_initial * alpha)) *
            wigner_uppercase_d_matrix<typename traits::NumTraits<Scalar>::real_t>(
                   f, m_initial, m_final, 0, beta, 0) *
            Scalar(std::cos(-m_final * gamma), std::sin(-m_final * gamma));
    } else {
        if (std::abs(
                std::remainder(m_initial * alpha, PI<typename traits::NumTraits<Scalar>::real_t>)) >
            numerical_precision) {
            throw std::invalid_argument(
                "The scalar type must be complex if m_initial*alpha is not a multiple of pi");
        }
        if (std::abs(
                std::remainder(m_final * gamma, PI<typename traits::NumTraits<Scalar>::real_t>)) >
            numerical_precision) {
            throw std::invalid_argument(
                "The scalar type must be complex if m_final*gamma is not a multiple of pi");
        }
        std::complex<Scalar> result = 0;
        for (typename traits::NumTraits<Scalar>::real_t m = -f;
             m <= f; // NOSONAR m is precisely representable
             ++m) {
            result += wigner_uppercase_d_matrix_pi_half(f, m_initial, m) *
                std::complex<typename traits::NumTraits<Scalar>::real_t>(std::cos(-m * beta),
                                                                         std::sin(-m * beta)) *
                wigner_uppercase_d_matrix_pi_half(f, m, -m_final);
        }
        result *= std::pow(std::complex<typename traits::NumTraits<Scalar>::real_t>(0, 1),
                           2 * f - m_initial - m_final) *
            static_cast<typename traits::NumTraits<Scalar>::real_t>(std::pow(-1, 2 * m_initial));
        return result.real();
    }
}

} // namespace pairinteraction::wigner
