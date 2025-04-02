// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/utils/wigner.hpp"

#include <complex>
#include <doctest/doctest.h>
#include <limits>

namespace pairinteraction {
DOCTEST_TEST_CASE("construction of wigner D matrix") {
    constexpr double PI = 3.141592653589793238462643383279502884;
    constexpr double numerical_precision = 100 * std::numeric_limits<double>::epsilon();

    auto wigner_real_entry =
        wigner::wigner_uppercase_d_matrix<double>(0.5, 0.5, -0.5, 4 * PI, PI / 3, 2 * PI);
    auto wigner_real_entry_reference = -0.5;
    DOCTEST_CHECK((wigner_real_entry - wigner_real_entry_reference) <= numerical_precision);

    std::string error_msg =
        "The scalar type must be complex if m_initial*alpha is not a multiple of pi";
    DOCTEST_CHECK_THROWS_WITH_AS(
        wigner::wigner_uppercase_d_matrix<double>(0.5, 0.5, -0.5, 0.1 * PI, 0, 0);
        , error_msg.c_str(), std::invalid_argument);

    error_msg = "The scalar type must be complex if m_final*gamma is not a multiple of pi";
    DOCTEST_CHECK_THROWS_WITH_AS(
        wigner::wigner_uppercase_d_matrix<double>(0.5, 0.5, -0.5, 0, 0, 0.1 * PI);
        , error_msg.c_str(), std::invalid_argument);

    auto wigner_complex_entry = wigner::wigner_uppercase_d_matrix<std::complex<double>>(
        0.5, 0.5, -0.5, 0.5 * PI, PI, -0.5 * PI);
    auto wigner_complex_entry_reference = std::complex<double>(0, 1);
    DOCTEST_CHECK(std::abs(wigner_complex_entry - wigner_complex_entry_reference) <=
                  numerical_precision);
}
} // namespace pairinteraction
