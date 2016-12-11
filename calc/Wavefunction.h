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

#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <string>
#include <vector>
#include "dtypes.h"
#include "QuantumDefect.h"

// --- Numerov's method ---

class Numerov {
    QuantumDefect const& qd;
    std::vector<real_t> x;
public:
    real_t const dx;
    Numerov(QuantumDefect const& qd);
    std::vector<real_t> axis() const;
    std::vector<real_t> integrate();

    constexpr static inline int power_kernel(int power)
    {
        return 2*power+2;
    }
};

// --- Whittaker method ---

class Whittaker {
    QuantumDefect const& qd;
    std::vector<real_t> x;
public:
    real_t const dx;
    Whittaker(QuantumDefect const& qd);
    std::vector<real_t> axis() const;
    std::vector<real_t> integrate();

    constexpr static inline real_t power_kernel(int power)
    {
        return 1.5*power;
    }
};


// --- Matrix element calculation ---


template < typename T >
size_t findidx(T const& x, typename T::value_type const& d) {
    auto it = std::find(std::begin(x), std::end(x), d);
    return std::distance(std::begin(x), it);
}


template < typename T >
real_t IntegrateRadialElement( QuantumDefect const& qd1, int power, QuantumDefect const& qd2)
{
    T N1(qd1);
    T N2(qd2);

    auto const& x1 = N1.axis();
    auto const& y1 = N1.integrate();
    auto const& x2 = N2.axis();
    auto const& y2 = N2.integrate();
    auto const  dx = N1.dx;

    auto const xmin = x1.front() >= x2.front() ? x1.front() : x2.front();
    auto const xmax = x1.back() <= x2.back() ? x1.back() : x2.back();

    real_t mu = 0;
    // If there is an overlap, calculate the matrix element
    if (xmin <= xmax) {
        int start1 = findidx(x1, xmin);
        int end1   = findidx(x1, xmax);
        int start2 = findidx(x2, xmin);
        int end2   = findidx(x2, xmax);

        int i1, i2;
        for (i1 = start1, i2 = start2; i1 < end1 && i2 < end2; ++i1, ++i2)
        {
            mu += y1[i1]*y2[i2] * std::pow(x1[i1], T::power_kernel(power)) * dx;
        }
        mu = std::abs(2*mu);
    }

    return mu;
}


#endif // WAVEFUNCTION_H
