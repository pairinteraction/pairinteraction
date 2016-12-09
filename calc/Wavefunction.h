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

#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

#include <string>
#include <vector>
#include "dtypes.h"
#include "QuantumDefect.h"

// --- Numerov's method ---

class Numerov {
  QuantumDefect const& qd;
  std::vector<real_t> x;
public:
  const real_t dx;
  Numerov(QuantumDefect const& qd);
  std::vector<real_t> axis();
  std::vector<real_t> integrate();
};

// --- Whittaker method ---

class Whittaker {
  QuantumDefect const& qd;
  std::vector<real_t> x;
public:
  const real_t dx;
  Whittaker(QuantumDefect const& qd);
  std::vector<real_t> axis();
  std::vector<real_t> integrate();
};

// --- Matrix element calculation ---

template < typename T >
size_t findidx(T x, real_t d) {
    auto it = std::find(std::begin(x), std::end(x), d);
    return std::distance(std::begin(x), it);
}

real_t IntegrateRadialElement(Numerov N1, int power, Numerov N2);

real_t IntegrateRadialElement(Whittaker N1, int power, Whittaker N2);

#endif // WAVEFUNCTION_HPP
