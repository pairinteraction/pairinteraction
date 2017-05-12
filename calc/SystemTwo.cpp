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
#include "SystemTwo.h"

#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>
#include <unordered_set>

SystemTwo::SystemTwo(const SystemOne &b1, const SystemOne &b2)
    : basis1(b1), basis2(b2)
{}

std::vector<StateOne> SystemTwo::getStatesFirst() {
    this->build();
    std::unordered_set<StateOne> states_one; // TODO make set work (this would have the benefit over unordered_set that the states are sorted)
    for (const auto &state : states) {
        states_one.insert(StateOne(state.getFirstState()));
    }
    return std::vector<StateOne>(states_one.begin(), states_one.end());
}

std::vector<StateOne> SystemTwo::getStatesSecond() {
    this->build();
    std::unordered_set<StateOne> states_one; // TODO make set work (this would have the benefit over unordered_set that the states are sorted)
    for (const auto &state : states) {
        states_one.insert(StateOne(state.getSecondState()));
    }
    return std::vector<StateOne>(states_one.begin(), states_one.end());
}

void SystemTwo::initialize()
{
    // Restrict one atom states to the allowed quantum numbers
    basis1.build();
    basis1.restrictN(range_n);
    basis1.restrictL(range_l);
    basis1.restrictJ(range_j);
    basis1.restrictM(range_m);
    basis2.build();
    basis2.restrictN(range_n);
    basis2.restrictL(range_l);
    basis2.restrictJ(range_j);
    basis2.restrictM(range_m);

    // TODO check basis1.getEnergies().size() == basis1.getCoefficients().outerSize() == basis1.getCoefficients().cols(), use method that returns the value of basis1.getEnergies().size()
    // TODO check basis1.getStates().size() == basis1.getCoefficients().innerSize() == basis1.getCoefficients().rows()
    // TODO check basis2.getEnergies().size() == basis2.getCoefficients().outerSize() == basis2.getCoefficients().cols()
    // TODO check basis2.getStates().size() == basis2.getCoefficients().innerSize() == basis2.getCoefficients().rows()

    // TODO consider symmetries

    // Combine the one atom states
    energies.reserve(basis1.getEnergies().size()*basis2.getEnergies().size());
    states.reserve(basis1.getStates().size()*basis2.getStates().size());
    std::vector<eigen_triplet_t> coefficients_triplets; // TODO reserve
    std::vector<double> sqnorm_list(basis1.getStates().size()*basis2.getStates().size(), 0);
    std::vector<size_t> stateidentifier2row(basis1.getStates().size()*basis2.getStates().size(), std::numeric_limits<size_t>::max());

    size_t col_new = 0;
    for (size_t col_1=0; col_1<basis1.getEnergies().size(); ++col_1) {
        for (size_t col_2=0; col_2<basis2.getEnergies().size(); ++col_2) {

            double energy = basis1.getEnergies()[col_1]+basis2.getEnergies()[col_2];
            if ((energy < energy_min && energy_min != std::numeric_limits<double_t>::lowest()) || (energy > energy_max  && energy_max != std::numeric_limits<double_t>::max())) continue;
            energies.push_back(energy);

            for (eigen_iterator_t triple_1(basis1.getCoefficients(),col_1); triple_1; ++triple_1) {
                for (eigen_iterator_t triple_2(basis2.getCoefficients(),col_2); triple_2; ++triple_2) {
                    size_t row_1 = triple_1.row();
                    size_t row_2 = triple_2.row();

                    scalar_t value_new = triple_1.value() * triple_2.value();

                    size_t stateidentifier = basis2.getStates().size()*row_1 + row_2;
                    size_t row_new = stateidentifier2row[stateidentifier];

                    if (row_new == std::numeric_limits<size_t>::max()) { // if stateidentifier not contained in map
                        row_new = states.size();
                        stateidentifier2row[stateidentifier] = row_new;
                        states.push_back(StateTwo(basis1.getStates()[row_1], basis2.getStates()[row_2]));
                    }

                    coefficients_triplets.push_back(eigen_triplet_t(row_new, col_new, value_new));
                    sqnorm_list[row_new] += std::pow(std::abs(value_new), 2);
                }
            }

            ++col_new;
        }
    }

    // Save storage
    basis1.clear();
    basis2.clear();

    states.shrink_to_fit();
    energies.shrink_to_fit();

    coefficients.resize(states.size(),energies.size());
    coefficients.setFromTriplets(coefficients_triplets.begin(), coefficients_triplets.end());
    coefficients_triplets.clear();

    // Build transformator and remove states (if the squared norm is to small) // TODO make this a method of the base class
    std::vector<StateTwo> states_new;
    states_new.reserve(states.size());
    std::vector<eigen_triplet_t> triplets_transformator;
    triplets_transformator.reserve(states.size());

    size_t idx_new = 0;
    for (size_t idx = 0; idx < states.size(); ++idx) {
        if (sqnorm_list[idx] > 0.05) {
            states_new.push_back(states[idx]);
            triplets_transformator.push_back(eigen_triplet_t(idx_new++,idx,1));
        }
    }

    states_new.shrink_to_fit();
    eigen_sparse_t transformator(idx_new,states.size());
    transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

    states = states_new;

    // Apply transformator in order to remove rows from the coefficient matrix (i.e. states)
    coefficients = transformator*coefficients;
}
