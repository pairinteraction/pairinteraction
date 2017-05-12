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

#include "dtypes.h"
#include "SystemOne.h"

#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>
#include <unordered_set>

SystemOne::SystemOne(std::string const& element)
    : element(element)
{}

const std::string& SystemOne::getElement() const {
    return element;
}

void SystemOne::initialize()
{
    // TODO check whether specified basis is finite

    // TODO consider symmetries

    size_t idx = 0;
    std::vector<eigen_triplet_t> coefficients_triplets; // TODO reserve states, energies, basis_triplets
    std::set<int> range_adapted_n, range_adapted_l;
    std::set<float> range_adapted_j, range_adapted_m;

    if (range_n.empty()) {
        range_adapted_n = std::set<int>({}); // TODO if empty, calculate the range via the energies
    } else {
        range_adapted_n = range_n;
    }
    for (auto n : range_adapted_n) {

        if (range_l.empty()) {
            for (int l = 0; l < n; ++l) {
                range_adapted_l.insert(l);
            }
        } else {
            range_adapted_l = range_l;
        }
        for (auto l : range_adapted_l) {
            if (l > n-1) continue;

            if (range_j.empty()) {
                range_adapted_j = std::set<float>({std::fabs(l-0.5f), l+0.5f});
            } else {
                range_adapted_j = range_j;
            }
            for (auto j : range_adapted_j) {
                if (std::fabs(j-l) != 0.5) continue;

                double energy = StateOne(element,n,l,j,0.5).getEnergy();
                if ((energy < energy_min && energy_min != std::numeric_limits<double>::lowest()) || (energy > energy_max  && energy_max != std::numeric_limits<double>::max())) continue; // TODO take into account numerical errors

                if (range_m.empty()) {
                    for (float m = -j; m < j; ++m) {
                        range_adapted_l.insert(m);
                    }
                } else {
                    range_adapted_m = range_m;
                }
                for (auto m : range_adapted_m) {
                    if (std::fabs(m) > j) continue;

                    states.push_back(StateOne(element,n,l,j,m));
                    energies.push_back(energy);
                    coefficients_triplets.push_back(eigen_triplet_t(idx,idx,1)); // TODO take into account symmetries

                    ++idx;
                }
            }
        }
    }

    coefficients.resize(idx,idx);
    coefficients.setFromTriplets(coefficients_triplets.begin(), coefficients_triplets.end());
    coefficients_triplets.clear();
}
