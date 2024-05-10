#pragma once

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

#include "utils/maths.hpp"
#include "utils/traits.hpp"

namespace wigner {

namespace {

template <typename Real>
inline constexpr Real PI = 3.141592653589793238462643383279502884;

template <typename Real>
inline Real wigner_uppercase_d_matrix_pi_half(Real f, Real m_initial, Real m_final) {
    Real result = 0;
    for (Real k = std::max(0, static_cast<int>(m_final - m_initial));
         k <= f + std::min(m_final, -m_initial); ++k) {
        result += std::pow(-1, k) * maths::binomial_coefficient(f + m_final, k) *
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
    if constexpr (traits::is_complex_v<Scalar>) {
        return Scalar(std::cos(-m_initial * alpha), std::sin(-m_initial * alpha)) *
            wigner_uppercase_d_matrix<typename traits::NumTraits<Scalar>::real_t>(
                   f, m_initial, m_final, 0, beta, 0) *
            Scalar(std::cos(-m_final * gamma), std::sin(-m_final * gamma));
    } else {
        if (std::abs(
                std::remainder(m_initial * alpha, PI<typename traits::NumTraits<Scalar>::real_t>)) >
            10 * std::numeric_limits<typename traits::NumTraits<Scalar>::real_t>::epsilon()) {
            throw std::invalid_argument(
                "The scalar type must be complex if m_initial*alpha is not a multiple of pi");
        }
        if (std::abs(
                std::remainder(m_final * gamma, PI<typename traits::NumTraits<Scalar>::real_t>)) >
            10 * std::numeric_limits<typename traits::NumTraits<Scalar>::real_t>::epsilon()) {
            throw std::invalid_argument(
                "The scalar type must be complex if m_final*gamma is not a multiple of pi");
        }
        std::complex<Scalar> result = 0;
        for (typename traits::NumTraits<Scalar>::real_t m = -f; m <= f; ++m) {
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

} // namespace wigner
