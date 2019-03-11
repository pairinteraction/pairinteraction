/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * This file is part of the pairinteraction library.
 *
 * The pairinteraction library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The pairinteraction library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "WignerD.h"

WignerD::WignerD() = default;

double WignerD::operator()(float j, float m, float mp, double beta) {
    double tolerance = 1e-16;
    if (std::abs(beta - M_PI / 2) < tolerance) {
        return this->evalWignerdPiHalf(j, m, mp);
    }
    return this->evalWignerd(j, m, mp, beta);
}

std::complex<double> WignerD::operator()(float j, float m, float mp, double alpha, double beta,
                                         double gamma) {
    return std::complex<double>(std::cos(-m * alpha), std::sin(-m * alpha)) *
        this->operator()(j, m, mp, beta) *
        std::complex<double>(std::cos(-mp * gamma), std::sin(-mp * gamma));
}

double WignerD::evalWignerdPiHalf(float j, float m, float mp) {
    double r = 0;
    for (unsigned int k = std::max(0, static_cast<int>(mp - m)); k <= j + std::min(mp, -m); ++k) {
        r += std::pow(-1, k) * boost::math::binomial_coefficient<double>(j + mp, k) *
            boost::math::binomial_coefficient<double>(j - mp, k + m - mp);
    }
    r *= std::pow(-1., m - mp) / std::pow(2., j) *
        std::sqrt(
             boost::math::factorial<double>(j + m) * boost::math::factorial<double>(j - m) /
             (boost::math::factorial<double>(j + mp) * boost::math::factorial<double>(j - mp)));
    return r;
}

double WignerD::evalWignerd(float j, float m, float mp, double beta) {
    std::complex<double> r = 0;
    for (float mpp = j; mpp >= -j; --mpp) { // TODO parallelize
        r += this->evalWignerdPiHalf(j, m, mpp) *
            std::complex<double>(std::cos(-mpp * beta), std::sin(-mpp * beta)) *
            this->evalWignerdPiHalf(j, mpp, -mp);
    }
    r *= std::pow(std::complex<double>(0, 1), 2. * j - m - mp) * std::pow(-1., 2. * m);
    return r.real();
}
