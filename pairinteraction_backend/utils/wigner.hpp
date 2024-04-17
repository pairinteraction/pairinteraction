#pragma once

#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <stdexcept>

#include "utils/maths.hpp"
#include "utils/traits.hpp"

namespace wigner {

namespace {

double wigner_uppercase_d_matrix_pi_half(float f, float m_initial, float m_final) {
    double result = 0;
    for (int k = std::max(0, static_cast<int>(m_final - m_initial));
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
requires(!traits::is_complex_v<Scalar>) Scalar
    wigner_uppercase_d_matrix(float f, float m_initial, float m_final, double alpha, double beta,
                              double gamma) {
    if (std::abs(std::remainder(m_initial * alpha, std::numbers::pi)) >
        10 * std::numeric_limits<double>::epsilon()) {
        throw std::invalid_argument(
            "The scalar type must be complex if m_initial*alpha is not a multiple of pi");
    }
    if (std::abs(std::remainder(m_final * gamma, std::numbers::pi)) >
        10 * std::numeric_limits<double>::epsilon()) {
        throw std::invalid_argument(
            "The scalar type must be complex if m_final*gamma is not a multiple of pi");
    }
    std::complex<Scalar> result = 0;
    for (float m = -f; m <= f; ++m) {
        result += wigner_uppercase_d_matrix_pi_half(f, m_initial, m) *
            std::complex<double>(std::cos(-m * beta), std::sin(-m * beta)) *
            wigner_uppercase_d_matrix_pi_half(f, m, -m_final);
    }
    result *= std::pow(std::complex<double>(0, 1), 2 * f - m_initial - m_final) *
        std::pow(-1, 2 * m_initial);
    return result.real();
}

template <typename Scalar>
requires traits::is_complex_v<Scalar> Scalar wigner_uppercase_d_matrix(float f, float m_initial,
                                                                       float m_final, double alpha,
                                                                       double beta, double gamma) {
    return Scalar(std::cos(-m_initial * alpha), std::sin(-m_initial * alpha)) *
        wigner_uppercase_d_matrix<typename traits::NumTraits<Scalar>::real_t>(f, m_initial, m_final,
                                                                              0, beta, 0) *
        Scalar(std::cos(-m_final * gamma), std::sin(-m_final * gamma));
}

} // namespace wigner
