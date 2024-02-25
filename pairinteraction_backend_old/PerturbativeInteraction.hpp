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

#ifndef PERTURBATIVEINTERACTION_H
#define PERTURBATIVEINTERACTION_H

#include "MatrixElementCache.hpp"
#include "State.hpp"

#include <array>
#include <vector>

class PerturbativeInteraction {
public:
    PerturbativeInteraction(MatrixElementCache &cache);
    PerturbativeInteraction(double angle, MatrixElementCache &cache);
    double getC6(const StateTwo &state, double deltaN); // return value in GHz*um^6
    Eigen::MatrixX<double> getC6(const std::vector<StateTwo> &states,
                                 double deltaN);                       // return value in GHz*um^6
    Eigen::MatrixX<double> getC3(const std::vector<StateTwo> &states); // return value in GHz*um^3
    Eigen::MatrixX<double> getEnergy(const std::vector<StateTwo> &states); // return value in GHz
private:
    void initializeAngleTerms(double angle);
    MatrixElementCache &cache;
    std::vector<std::array<int, 2>> array_q;
    std::array<double, 9> array_angle_term;
};

#endif
