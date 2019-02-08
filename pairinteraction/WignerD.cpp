/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
