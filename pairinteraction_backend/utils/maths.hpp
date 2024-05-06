#pragma once

#include <stdexcept>

namespace maths {

inline double binomial_coefficient(double n, double k) {
    if (n < k || k < 0) {
        throw std::invalid_argument("It must be n >= k >= 0.");
    }
    if (k == 0) {
        return 1;
    }

    // Recursion (n, k) = (n - 1, k - 1) * n / k
    double result = n - k + 1;
    for (int i = 1; i < k; ++i) {
        result *= (n - k + 1 + i) / (i + 1);
    }

    return result;
}

inline double factorial(double n) {
    if (n < 0) {
        throw std::invalid_argument("It must be n >= 0.");
    }
    if (n == 0) {
        return 1;
    }

    double result = 1;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }

    return result;
}

} // namespace maths
