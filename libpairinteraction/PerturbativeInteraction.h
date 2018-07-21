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

#ifndef PERTURBATIVEINTERACTION_H
#define PERTURBATIVEINTERACTION_H

#include "MatrixElementCache.h"
#include "State.h"
#include "dtypes.h"

#include <array>
#include <vector>

class PerturbativeInteraction {
public:
    PerturbativeInteraction(MatrixElementCache &cache);
    PerturbativeInteraction(double angle, MatrixElementCache &cache);
    PerturbativeInteraction(double angle, double weak_bfield_along_z, MatrixElementCache &cache);
    double getC6(StateTwo state, double deltaN); // return value in GHz*um^6
    eigen_dense_double_t getC6(std::vector<StateTwo> states,
                               double deltaN);                      // return value in GHz*um^6
    eigen_dense_double_t getC3(std::vector<StateTwo> states);       // return value in GHz*um^3
    eigen_dense_double_t getEnergies(std::vector<StateTwo> states); // return value in GHz
private:
    void initializeAngleTerms(double angle);
    MatrixElementCache &cache;
    double bfield;
    std::vector<std::array<int, 2>> array_q;
    std::array<double, 9> array_angle_term;
};

#endif
