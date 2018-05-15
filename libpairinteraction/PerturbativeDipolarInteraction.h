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

#ifndef PERTURBATIVEDIPOLARINTERACTION_H
#define PERTURBATIVEDIPOLARINTERACTION_H

#include "dtypes.h"
#include "MatrixElementCache.h"
#include "State.h"

#include <vector>
#include <array>

class PerturbativeDipolarInteraction {
public:
    PerturbativeDipolarInteraction(MatrixElementCache &cache);
    PerturbativeDipolarInteraction(MatrixElementCache &cache, double angle);
    double getC6(StateTwo state, double deltaN); // return value in GHz*um^6
    eigen_dense_double_t getC6(std::vector<StateTwo> states, double deltaN); // return value in GHz*um^6
    eigen_dense_double_t getC3(std::vector<StateTwo> states); // return value in GHz*um^3
private:
    MatrixElementCache &cache;
    std::vector<std::array<int,2>> array_q;
    std::array<double, 9> array_angle_term;
};

#endif
