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

SystemTwo::SystemTwo(const SystemOne &b1, const SystemOne &b2, std::string cachedir)
    : SystemBase(cachedir), element({{b1.getElement(), b2.getElement()}}), system1(b1), system2(b2), distance(std::numeric_limits<double>::max()), angle(0), kappa_max(1), sym_permutation(NA) {
}

SystemTwo::SystemTwo(const SystemOne &b1, const SystemOne &b2, std::string cachedir, bool memory_saving)
    : SystemBase(cachedir, memory_saving), element({{b1.getElement(), b2.getElement()}}), system1(b1), system2(b2), distance(std::numeric_limits<double>::max()), angle(0), kappa_max(1), sym_permutation(NA) {
}

SystemTwo::SystemTwo(const SystemOne &b1, const SystemOne &b2)
    : SystemBase(), element({{b1.getElement(), b2.getElement()}}), system1(b1), system2(b2), distance(std::numeric_limits<double>::max()), angle(0), kappa_max(1), sym_permutation(NA) {
}

SystemTwo::SystemTwo(const SystemOne &b1, const SystemOne &b2, bool memory_saving)
    : SystemBase(memory_saving), element({{b1.getElement(), b2.getElement()}}), system1(b1), system2(b2), distance(std::numeric_limits<double>::max()), angle(0), kappa_max(1), sym_permutation(NA) {
}

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

const std::array<std::string, 2>& SystemTwo::getElement() {
    return element;
}

void SystemTwo::setDistance(double d) {
    this->onParameterChange();
    distance = d;
}

void SystemTwo::setAngle(double a) {
    this->onParameterChange();
    angle = a;
}

void SystemTwo::setConservedParityUnderPermutation(parity_t parity) {
    this->onSymmetryChange();
    sym_permutation = parity;
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to initialize Basis //////////////
////////////////////////////////////////////////////////////////////

void SystemTwo::initializeBasis()
{
    ////////////////////////////////////////////////////////////////////
    /// Restrict one atom states to the allowed quantum numbers ////////
    ////////////////////////////////////////////////////////////////////

    system1.diagonalize(); // important!
    system1.restrictN(range_n);
    system1.restrictL(range_l);
    system1.restrictJ(range_j);
    system1.restrictM(range_m);
    system2.diagonalize(); // important!
    system2.restrictN(range_n);
    system2.restrictL(range_l);
    system2.restrictJ(range_j);
    system2.restrictM(range_m);

    ////////////////////////////////////////////////////////////////////
    /// Check wther the single atom states fit to the symmetries ///////
    ////////////////////////////////////////////////////////////////////

    if (sym_permutation != NA) {
        // TODO check system1 == system2
    }

    ////////////////////////////////////////////////////////////////////
    /// Combine one atom states ////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // TODO consider further symmetries and check whether they are applicable

    std::vector<eigen_triplet_t>  hamiltonianmatrix_triplets;
    hamiltonianmatrix_triplets.reserve(system1.getNumVectors()*system2.getNumVectors());
    states.reserve(system1.getNumStates()*system2.getNumStates());
    std::vector<eigen_triplet_t> coefficients_triplets; // TODO reserve
    std::vector<double> sqnorm_list(system1.getNumStates()*system2.getNumStates(), 0);
    std::vector<size_t> stateidentifier2row(system1.getNumStates()*system2.getNumStates(), std::numeric_limits<size_t>::max());

    size_t col_new = 0;
    for (size_t col_1=0; col_1<system1.getNumVectors(); ++col_1) {
        for (size_t col_2=0; col_2<system2.getNumVectors(); ++col_2) {

            // In case of permutation symmetry: skip half of the basis vector pairs
            if ((sym_permutation == EVEN && col_1 <= col_2) || (sym_permutation == ODD && col_1 < col_2)) { // asym
                continue;
            }

            // Continue if the pair statet energy is not valid
            double energy = this->real(system1.getHamiltonianmatrix().coeff(col_1, col_1)+system2.getHamiltonianmatrix().coeff(col_2, col_2));
            if (!checkIsEnergyValid(energy)) continue;

            // Store the pair state energy
            hamiltonianmatrix_triplets.push_back(eigen_triplet_t(col_new, col_new, energy));

            // Build the basis vector that corresponds to the stored pair state energy
            for (eigen_iterator_t triple_1(system1.getCoefficients(),col_1); triple_1; ++triple_1) {
                for (eigen_iterator_t triple_2(system2.getCoefficients(),col_2); triple_2; ++triple_2) {
                    size_t row_1 = triple_1.row();
                    size_t row_2 = triple_2.row();

                    scalar_t value_new = triple_1.value() * triple_2.value();

                    // Adapt the normalization if required if required by symmetries
                    if (sym_permutation != NA && col_1 != col_2 ) {
                        value_new /= std::sqrt(2);
                    }

                    // Add an entry to the current basis vector
                    this->addCoefficient(row_1, row_2, col_new, value_new, stateidentifier2row, coefficients_triplets, sqnorm_list);

                    // Add further entries to the current basis vector if required by symmetries
                    if (sym_permutation != NA && col_1 != col_2 ) {
                        value_new *= (sym_permutation == EVEN) ? -1 : 1;
                        this->addCoefficient(row_2, row_1, col_new, value_new, stateidentifier2row, coefficients_triplets, sqnorm_list);
                    }
                }
            }

            ++col_new;
        }
    }

    // Delete unecessary storage
    system1 = SystemOne(element[0], cachedir.string());
    system2 = SystemOne(element[1], cachedir.string());

    // Build data
    states.shrink_to_fit();

    coefficients.resize(states.size(),col_new);
    coefficients.setFromTriplets(coefficients_triplets.begin(), coefficients_triplets.end());
    coefficients_triplets.clear();

    hamiltonianmatrix.resize(col_new,col_new);
    hamiltonianmatrix.setFromTriplets(hamiltonianmatrix_triplets.begin(), hamiltonianmatrix_triplets.end());
    hamiltonianmatrix_triplets.clear();

    ////////////////////////////////////////////////////////////////////
    /// Remove states that barely occur ////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // Build transformator and remove states if the squared norm is to small
    removeRestrictedStates([=](size_t idx) -> bool { return sqnorm_list[idx] > 0.05; } );
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to calculate the interaction /////
////////////////////////////////////////////////////////////////////

void SystemTwo::initializeInteraction() {
    if (distance == std::numeric_limits<double>::max()) {
        interaction.clear();
        return;
    }

    ////////////////////////////////////////////////////////////////////
    /// Prepare the calculation of the interaction /////////////////////
    ////////////////////////////////////////////////////////////////////

    // Check if something to do
    std::vector<bool> calculation_required(4, false); // TODO use array or set

    double tolerance = 1e-24;

    if (interaction.find(0) == interaction.end()) {
        calculation_required[0] = true;
    }

    if (std::fabs(1-3*std::pow(std::cos(angle),2)) > tolerance && interaction.find(1) == interaction.end()) {
        calculation_required[1] = true;
    }

    if (std::fabs(1.5*std::pow(std::sin(angle),2)) > tolerance && interaction.find(2) == interaction.end()) {
        calculation_required[2] = true;
    }

    if (std::fabs(3./std::sqrt(2)*std::sin(angle)*std::cos(angle)) > tolerance && interaction.find(3) == interaction.end()) {
        calculation_required[3] = true;
    }

    if (!std::accumulate(calculation_required.begin(), calculation_required.end(), false)) return;

    // TODO add operators for higer order interaction !!!

    // Precalculate matrix elements
    std::string matrixelementsdir = "";
    if (!cachedir.empty()) matrixelementsdir = (cachedir / "cache_elements.db").string(); // TODO do this in the MatrixElements class, just pass cachedir as an argument to the constructor

    MatrixElements matrixelements1(element[0], matrixelementsdir);
    MatrixElements matrixelements2(element[1], matrixelementsdir);

    auto states1 = this->getStatesFirst();
    auto states2 = this->getStatesSecond();

    for (int kappa = 1; kappa <= kappa_max; ++kappa) {
        matrixelements1.precalculateMultipole(states1, kappa);
        matrixelements2.precalculateMultipole(states2, kappa);
    }

    ////////////////////////////////////////////////////////////////////
    /// Generate the interaction in the canonical basis ////////////////
    ////////////////////////////////////////////////////////////////////

    std::unordered_map<int, std::vector<eigen_triplet_t>> interaction_triplets; // TODO reserve

    for (size_t col=0; col<states.size(); ++col) { // TODO parallelization
        for (size_t row=0; row<states.size(); ++row) {
            //if (row < col) continue; // TODO use this restriction and construct the total matrix (that is needed in order to be transformable) from the diagonal matrix afterwards

            const StateTwo &state_row = states[row];
            const StateTwo &state_col = states[col];

            if (state_row.element.empty() || state_col.element.empty()  ) continue; // TODO artifical states TODO [dummystates]

            if (selectionRulesMultipole(state_row.first(), state_col.first(), 1) && selectionRulesMultipole(state_row.second(), state_col.second(), 1)) {
                int q1 = state_row.first().m-state_col.first().m;
                int q2 = state_row.second().m-state_col.second().m;

                if (q1 == 0 && q2 == 0 && calculation_required[1]) {
                    scalar_t val = inverse_electric_constant * matrixelements1.getMultipole(state_row.first(), state_col.first(), 1) *
                            matrixelements2.getMultipole(state_row.second(), state_col.second(), 1);
                    interaction_triplets[1].push_back(eigen_triplet_t(row, col, val));
                } else if (q1 != 0 && q2 != 0 && q1+q2 == 0 && (calculation_required[0] || calculation_required[2])) {
                    scalar_t val = inverse_electric_constant * matrixelements1.getMultipole(state_row.first(), state_col.first(), 1) *
                            matrixelements2.getMultipole(state_row.second(), state_col.second(), 1);
                    if (calculation_required[0]) interaction_triplets[0].push_back(eigen_triplet_t(row, col, val));
                    if (calculation_required[2]) interaction_triplets[2].push_back(eigen_triplet_t(row, col, -val));
                } else if (std::abs(q1+q2) == 1 && calculation_required[3]) {
                    scalar_t val = inverse_electric_constant * matrixelements1.getMultipole(state_row.first(), state_col.first(), 1) *
                            matrixelements2.getMultipole(state_row.second(), state_col.second(), 1);
                    if (q1 == 1 || q2 == 1) val *= -1;// TODO think of a better way
                    interaction_triplets[3].push_back(eigen_triplet_t(row, col, val));
                } else if (std::abs(q1+q2) == 2 && calculation_required[2]) {
                    scalar_t val = inverse_electric_constant * matrixelements1.getMultipole(state_row.first(), state_col.first(), 1) *
                            matrixelements2.getMultipole(state_row.second(), state_col.second(), 1);
                    interaction_triplets[2].push_back(eigen_triplet_t(row, col, val));
                } // TODO simplify code

                // TODO add operators for higer order interaction !!!
            }
        }
    }

    for (size_t i = 0; i < calculation_required.size(); ++i) {
        if (!calculation_required[i]) continue;
        interaction[i].resize(states.size(),states.size());
        interaction[i].setFromTriplets(interaction_triplets[i].begin(), interaction_triplets[i].end());
        interaction_triplets[i].clear();
    }

    ////////////////////////////////////////////////////////////////////
    /// Transform the interaction to the used basis ////////////////////
    ////////////////////////////////////////////////////////////////////

    for (auto &entry : interaction) entry.second = coefficients.adjoint()*entry.second*coefficients;
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to construct Hamiltonian /////////
////////////////////////////////////////////////////////////////////

void SystemTwo::addInteraction() {
    if (distance == std::numeric_limits<double>::max()) return;

    // Calculate the distance dependency
    double powerlaw = 1./std::pow(distance, 3);

    // Build the total Hamiltonian
    double tolerance = 1e-24; // TODO think about it, can the expressions be evalualted identical to zero?
    if (interaction.find(0) != interaction.end()) hamiltonianmatrix -= powerlaw*interaction[0];
    if (std::fabs(1-3*std::pow(std::cos(angle),2)) > tolerance) hamiltonianmatrix += (1-3*std::pow(std::cos(angle),2))*powerlaw*interaction[1];
    if (std::fabs(1.5*std::pow(std::sin(angle),2)) > tolerance) hamiltonianmatrix -= (1.5*std::pow(std::sin(angle),2))*powerlaw*interaction[2];
    if (std::fabs(3./std::sqrt(2)*std::sin(angle)*std::cos(angle)) > tolerance) hamiltonianmatrix -= (3./std::sqrt(2)*std::sin(angle)*std::cos(angle))*powerlaw*interaction[3];

    // TODO add operators for higer order interaction !!!
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to transform the interaction /////
////////////////////////////////////////////////////////////////////

void SystemTwo::transformInteraction(const eigen_sparse_t &transformator)  {
    for (auto &entry : interaction) entry.second = transformator.adjoint()*entry.second*transformator;
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to delete the interaction ////////
////////////////////////////////////////////////////////////////////

void SystemTwo::deleteInteraction()  {
    interaction.clear();
}

////////////////////////////////////////////////////////////////////
/// Utility methods ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

void SystemTwo::addCoefficient(const size_t &row_1, const size_t &row_2, const size_t &col_new, const scalar_t &value_new, std::vector<size_t> &stateidentifier2row, std::vector<eigen_triplet_t> &coefficients_triplets, std::vector<double> &sqnorm_list) {
    size_t stateidentifier = system2.getNumStates()*row_1 + row_2;
    size_t row_new = stateidentifier2row[stateidentifier];

    if (row_new == std::numeric_limits<size_t>::max()) { // if stateidentifier not contained in map
        row_new = states.size();
        stateidentifier2row[stateidentifier] = row_new;
        states.push_back(StateTwo(system1.getStates()[row_1], system2.getStates()[row_2]));
    }

    coefficients_triplets.push_back(eigen_triplet_t(row_new, col_new, value_new));
    sqnorm_list[row_new] += std::pow(std::abs(value_new), 2);
}
