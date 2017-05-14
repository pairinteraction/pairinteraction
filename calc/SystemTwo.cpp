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
#include "MatrixElements.h"

#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>
#include <unordered_set>
#include <complex>

SystemTwo::SystemTwo(const SystemOne &b1, const SystemOne &b2, std::string cachedir)
    : SystemBase(cachedir), basis1(b1), basis2(b2)
{}

std::vector<StateOne> SystemTwo::getStatesFirst() {
    this->buildBasis();
    std::unordered_set<StateOne> states_one; // TODO make set work (this would have the benefit over unordered_set that the states are sorted)
    for (const auto &state : states) {
        states_one.insert(StateOne(state.getFirstState()));
    }
    return std::vector<StateOne>(states_one.begin(), states_one.end());
}

std::vector<StateOne> SystemTwo::getStatesSecond() {
    this->buildBasis();
    std::unordered_set<StateOne> states_one; // TODO make set work (this would have the benefit over unordered_set that the states are sorted)
    for (const auto &state : states) {
        states_one.insert(StateOne(state.getSecondState()));
    }
    return std::vector<StateOne>(states_one.begin(), states_one.end());
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to initialize Basis //////////////
////////////////////////////////////////////////////////////////////

void SystemTwo::initializeBasis()
{
    ////////////////////////////////////////////////////////////////////
    /// Restrict one atom states to the allowed quantum numbers ////////
    ////////////////////////////////////////////////////////////////////

    basis1.diagonalize(); // important!
    basis1.restrictN(range_n);
    basis1.restrictL(range_l);
    basis1.restrictJ(range_j);
    basis1.restrictM(range_m);
    basis2.diagonalize(); // important!
    basis2.restrictN(range_n);
    basis2.restrictL(range_l);
    basis2.restrictJ(range_j);
    basis2.restrictM(range_m);

    ////////////////////////////////////////////////////////////////////
    /// Combine one atom states ////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // TODO consider symmetries

    std::vector<eigen_triplet_t>  hamiltonianmatrix_triplets;
    hamiltonianmatrix_triplets.reserve(basis1.getNumVectors()*basis2.getNumVectors());
    states.reserve(basis1.getNumStates()*basis2.getNumStates());
    std::vector<eigen_triplet_t> coefficients_triplets; // TODO reserve
    std::vector<double> sqnorm_list(basis1.getNumStates()*basis2.getNumStates(), 0);
    std::vector<size_t> stateidentifier2row(basis1.getNumStates()*basis2.getNumStates(), std::numeric_limits<size_t>::max());

    size_t col_new = 0;
    for (size_t col_1=0; col_1<basis1.getNumVectors(); ++col_1) {
        for (size_t col_2=0; col_2<basis2.getNumVectors(); ++col_2) {

            double energy = std::complex<double>(basis1.getHamiltonianmatrix().coeff(col_1, col_1)+basis2.getHamiltonianmatrix().coeff(col_2, col_2)).real();
            if (!checkIsEnergyValid(energy)) continue;
            hamiltonianmatrix_triplets.push_back(eigen_triplet_t(col_new, col_new, energy));

            for (eigen_iterator_t triple_1(basis1.getCoefficients(),col_1); triple_1; ++triple_1) {
                for (eigen_iterator_t triple_2(basis2.getCoefficients(),col_2); triple_2; ++triple_2) {
                    size_t row_1 = triple_1.row();
                    size_t row_2 = triple_2.row();

                    scalar_t value_new = triple_1.value() * triple_2.value();

                    size_t stateidentifier = basis2.getNumStates()*row_1 + row_2;
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

    // Delete unecessary storage
    basis1.clear();
    basis2.clear();

    // Build data
    states.shrink_to_fit();

    coefficients.resize(states.size(),col_new);
    coefficients.setFromTriplets(coefficients_triplets.begin(), coefficients_triplets.end());
    coefficients_triplets.clear();

    hamiltonianmatrix.resize(states.size(),col_new);
    hamiltonianmatrix.setFromTriplets(hamiltonianmatrix_triplets.begin(), hamiltonianmatrix_triplets.end());
    hamiltonianmatrix_triplets.clear();

    ////////////////////////////////////////////////////////////////////
    /// Remove states that barely occur ////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // Build transformator and remove states if the squared norm is to small
    removeRestrictedStates([=](size_t idx) -> bool { return sqnorm_list[idx] > 0.05; } );
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to initialize Hamiltonian helpers
////////////////////////////////////////////////////////////////////

void SystemTwo::initializeHamiltonianhelpers() {
    // TODO
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to construct Hamiltonian /////////
////////////////////////////////////////////////////////////////////

void SystemTwo::initializeHamiltonian() {
    // TODO
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to transform Hamiltonian helpers /
////////////////////////////////////////////////////////////////////

void SystemTwo::transformHamiltonianhelpers(const eigen_sparse_t &transformator)  {
    (void) transformator;
    // TODO
}

