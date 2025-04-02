// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <stdexcept>
#include <type_traits>

namespace pairinteraction::maths {

template <typename Real>
inline Real binomial_coefficient(Real n, Real k) {
    static_assert(std::is_arithmetic_v<Real>);

    if (n < k || k < 0) {
        throw std::invalid_argument("It must be n >= k >= 0.");
    }
    if (k == 0) {
        return 1;
    }

    // Recursion (n, k) = (n - 1, k - 1) * n / k
    Real result = n - k + 1;
    for (int i = 1; i < k; ++i) {
        result *= (n - k + 1 + i) / (i + 1);
    }

    return result;
}

template <typename Real>
inline Real factorial(Real n) {
    static_assert(std::is_arithmetic_v<Real>);

    if (n < 0) {
        throw std::invalid_argument("It must be n >= 0.");
    }
    if (n == 0) {
        return 1;
    }

    Real result = 1;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }

    return result;
}

} // namespace pairinteraction::maths
