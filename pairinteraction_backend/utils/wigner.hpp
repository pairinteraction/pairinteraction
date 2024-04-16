#pragma once

#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <stdexcept>

#include "utils/maths.hpp"

namespace wigner {

namespace {

template <typename real_t>
real_t wigner_uppercase_d_matrix_pi_half(real_t f, real_t m_initial, real_t m_final) {
    real_t result = 0;
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

template <typename real_t>
real_t wigner_uppercase_d_matrix(real_t f, real_t m_initial, real_t m_final, real_t alpha,
                                 real_t beta, real_t gamma) {
    if (std::abs(std::remainder(m_initial * alpha, std::numbers::pi)) >
        10 * std::numeric_limits<real_t>::epsilon()) {
        throw std::invalid_argument(
            "The scalar type must be complex if m_initial*alpha is not a multiple of pi");
    }
    if (std::abs(std::remainder(m_final * gamma, std::numbers::pi)) >
        10 * std::numeric_limits<real_t>::epsilon()) {
        throw std::invalid_argument(
            "The scalar type must be complex if m_final*gamma is not a multiple of pi");
    }
    std::complex<real_t> result = 0;
    for (real_t m = -f; m <= f; ++m) {
        result += wigner_uppercase_d_matrix_pi_half(f, m_initial, m) *
            std::complex<real_t>(std::cos(-m * beta), std::sin(-m * beta)) *
            wigner_uppercase_d_matrix_pi_half(f, m, -m_final);
    }
    result *= std::pow(std::complex<real_t>(0, 1), 2 * f - m_initial - m_final) *
        std::pow(-1, 2 * m_initial);
    return result.real();
}

template <typename real_t>
std::complex<real_t> wigner_uppercase_d_matrix(real_t f, real_t m_initial, real_t m_final,
                                               real_t alpha, real_t beta, real_t gamma) {
    return std::complex<real_t>(std::cos(-m_initial * alpha), std::sin(-m_initial * alpha)) *
        wigner_uppercase_d_matrix<real_t>(f, m_initial, m_final, 0, beta, 0) *
        std::complex<real_t>(std::cos(-m_final * gamma), std::sin(-m_final * gamma));
}

} // namespace wigner
