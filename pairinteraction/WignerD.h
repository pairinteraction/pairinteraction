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

#ifndef WIGNERD_H
#define WIGNERD_H

#define _USE_MATH_DEFINES

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <cmath>
#include <complex>
#include <numeric>
#include <vector>

class WignerD {
public:
    WignerD();
    double operator()(float j, float m, float mp, double beta);
    std::complex<double> operator()(float j, float m, float mp, double alpha, double beta,
                                    double gamma);

private:
    double evalWignerdPiHalf(float j, float m, float mp);
    double evalWignerd(float j, float m, float mp, double beta);
};

#endif
