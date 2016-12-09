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

#ifndef QUANTUM_DEFECT_HPP
#define QUANTUM_DEFECT_HPP

#include <string>
#include "dtypes.h"

//typedef double real_t;

class QuantumDefect {
private:
    std::string species_;
    int n_, l_;
    real_t j_;
    real_t ac_;
    int Z_;
    real_t a1_, a2_, a3_, a4_;
    real_t rc_;
    real_t nstar_;
    real_t energy_;
public:
    QuantumDefect(std::string species, int n, int l, real_t j);
    const std::string &species;
    const int &n, &l;
    const real_t &j;
    const real_t &ac;
    const int &Z;
    const real_t &a1, &a2, &a3, &a4;
    const real_t &rc;
    const real_t &nstar;
    const real_t &energy;
};

real_t energy_level(std::string species, int n, int l, real_t j);

#endif // QUANTUM_DEFECT_HPP
