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

SystemTwo::SystemTwo(const SystemOne &b1, const SystemOne &b2, std::wstring cachedir)
    : SystemBase(cachedir), element({{b1.getElement(), b2.getElement()}}), system1(b1), system2(b2), distance(std::numeric_limits<double>::max()), angle(0), ordermax(3), sym_permutation(NA) {
}

SystemTwo::SystemTwo(const SystemOne &b1, const SystemOne &b2, std::wstring cachedir, bool memory_saving)
    : SystemBase(cachedir, memory_saving), element({{b1.getElement(), b2.getElement()}}), system1(b1), system2(b2), distance(std::numeric_limits<double>::max()), angle(0), ordermax(3), sym_permutation(NA) {
}

SystemTwo::SystemTwo(const SystemOne &b1, const SystemOne &b2)
    : SystemBase(), element({{b1.getElement(), b2.getElement()}}), system1(b1), system2(b2), distance(std::numeric_limits<double>::max()), angle(0), ordermax(3), sym_permutation(NA) {
}

SystemTwo::SystemTwo(const SystemOne &b1, const SystemOne &b2, bool memory_saving)
    : SystemBase(memory_saving), element({{b1.getElement(), b2.getElement()}}), system1(b1), system2(b2), distance(std::numeric_limits<double>::max()), angle(0), ordermax(3), sym_permutation(NA) {
}

std::vector<StateOne> SystemTwo::getStatesFirst() {  // TODO @hmenke typemap for "state_set<StateOne>"
    this->buildBasis();
    std::unordered_set<StateOne> states_one;
    for (const auto &entry : states) {
        states_one.insert(StateOne(entry.state.getFirstState()));
    }
    return std::vector<StateOne>(states_one.begin(), states_one.end());
}

std::vector<StateOne> SystemTwo::getStatesSecond() { // TODO @hmenke typemap for "state_set<StateOne>"
    this->buildBasis();
    std::unordered_set<StateOne> states_one;
    for (const auto &entry : states) {
        states_one.insert(StateOne(entry.state.getSecondState()));
    }
    return std::vector<StateOne>(states_one.begin(), states_one.end());
}

const std::array<std::wstring, 2>& SystemTwo::getElement() {
    return element;
}

void SystemTwo::setDistance(double d) {
    this->onParameterChange();
    distance = d;
}

void SystemTwo::setAngle(double a) {
    if (a != 0 && ordermax > 3) throw std::runtime_error( "A non-zero interaction angle can be directly used only for dipole-dipole interaction.");

    this->onParameterChange();
    angle = a;

    angle_terms[0] = -1.;
    angle_terms[1] = 1.-3.*std::pow(std::cos(angle),2);
    angle_terms[2] = -1.5*std::pow(std::sin(angle),2);
    angle_terms[3] = -3./std::sqrt(2)*std::sin(angle)*std::cos(angle);
}

void SystemTwo::setOrder(double o) {
    if (angle != 0 && o > 3) throw std::runtime_error( "A non-zero interaction angle can be directly used only for dipole-dipole interaction.");

    this->onParameterChange();
    ordermax = o;
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

    system1.diagonalize(); // it is important to call this method here!
    system1.restrictN(range_n);
    system1.restrictL(range_l);
    system1.restrictJ(range_j);
    system1.restrictM(range_m);
    system2.diagonalize(); // it is important to call this method here!
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
                    this->addCoefficient(row_1, row_2, col_new, value_new, coefficients_triplets, sqnorm_list);

                    // Add further entries to the current basis vector if required by symmetries
                    if (sym_permutation != NA && col_1 != col_2 ) {
                        value_new *= (sym_permutation == EVEN) ? -1 : 1;
                        this->addCoefficient(row_2, row_1, col_new, value_new, coefficients_triplets, sqnorm_list);
                    }
                }
            }

            ++col_new;
        }
    }

    // Delete unecessary storage
    system1 = SystemOne(element[0]);
    system2 = SystemOne(element[1]);

    // Build data
    states.shrink_to_fit();

    coefficients.resize(states.size(),col_new);
    coefficients.setFromTriplets(coefficients_triplets.begin(), coefficients_triplets.end());
    coefficients_triplets.clear();

    hamiltonianmatrix.resize(col_new,col_new);
    hamiltonianmatrix.setFromTriplets(hamiltonianmatrix_triplets.begin(), hamiltonianmatrix_triplets.end());
    hamiltonianmatrix_triplets.clear();

    ////////////////////////////////////////////////////////////////////
    /// Remove vectors with too small norm /////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // Build transformator and remove vectors (if the squared norm is too small)
    std::vector<eigen_triplet_t> triplets_transformator;
    triplets_transformator.reserve(coefficients.cols());

    size_t idx_new = 0;
    for (int idx=0; idx<coefficients.cols(); ++idx) { // idx = col = num basis vector
        double_t sqnorm = 0;

        // Calculate the square norm of the columns of the coefficient matrix
        for (eigen_iterator_t triple(coefficients,idx); triple; ++triple) {
            sqnorm += std::pow(std::abs(triple.value()),2);
        }
        if (sqnorm > 0.05) {
            triplets_transformator.push_back(eigen_triplet_t(idx,idx_new++,1));
        }
    }

    this->applyRightsideTransformator(triplets_transformator);

    ////////////////////////////////////////////////////////////////////
    /// Remove states that barely occur ////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // Build transformator and remove states if the squared norm is to small
    removeRestrictedStates([=](const enumerated_state<StateTwo> &entry) -> bool { return sqnorm_list[entry.idx] > 0.05; } );

}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to calculate the interaction /////
////////////////////////////////////////////////////////////////////

void SystemTwo::initializeInteraction() {
    if (distance == std::numeric_limits<double>::max()) return;

    ////////////////////////////////////////////////////////////////////
    /// Prepare the calculation of the interaction /////////////////////
    ////////////////////////////////////////////////////////////////////

    // Check if something to do
    double tolerance = 1e-24;

    std::vector<bool> calculation_required(4, false);
    std::vector<int> orange;

    if (angle != 0) { // setAngle and setOrder take care that a non-zero angle cannot occur for other interaction than dipole-dipole

        for (size_t i = 0; i < 4; ++i) {
            if (std::abs(angle_terms[i]) > tolerance && interaction_angulardipole.find(i) == interaction_angulardipole.end()) {
                calculation_required[i] = true;
            }
        }

    } else {

        for (unsigned int order = 3; order <= ordermax; ++order) {
            if (interaction_multipole.find(order) == interaction_multipole.end()) {
                orange.push_back(order);
            }
        }

    }

    // Precalculate matrix elements
    std::string matrixelementsdir = "";
    if (!cachedir.empty()) matrixelementsdir = (cachedir / "cache_elements.db").string(); // TODO do this in the MatrixElements class, just pass cachedir as an argument to the constructor

    std::string tmp(element[0].begin(), element[0].end()); // TODO think of a better solution
    MatrixElements matrixelements1(tmp, matrixelementsdir);
    tmp = std::string(element[1].begin(), element[1].end());
    MatrixElements matrixelements2(tmp, matrixelementsdir);

    auto states1 = this->getStatesFirst();
    auto states2 = this->getStatesSecond();

    for (unsigned int kappa = 1; kappa <= ordermax-2; ++kappa) {
        matrixelements1.precalculateMultipole(states1, kappa);
        matrixelements2.precalculateMultipole(states2, kappa); // TODO check whether system1 == system2
    }

    ////////////////////////////////////////////////////////////////////
    /// Generate the interaction in the canonical basis ////////////////
    ////////////////////////////////////////////////////////////////////

    std::unordered_map<int, std::vector<eigen_triplet_t>> interaction_angulardipole_triplets; // TODO reserve
    std::unordered_map<int, std::vector<eigen_triplet_t>> interaction_multipole_triplets; // TODO reserve

    // Loop over column entries
    for (const auto &c: states) { // TODO parallelization

        if (c.state.element.empty()) continue; // TODO artifical states TODO [dummystates]

        // Loop over row entries
        for (const auto &r: states) {

            if (r.state.element.empty()) continue; // TODO artifical states TODO [dummystates]
            if (r.idx < c.idx) continue;

            int q1 = r.state.first().m-c.state.first().m;
            int q2 = r.state.second().m-c.state.second().m;

            if (angle != 0) { // setAngle and setOrder take care that a non-zero angle cannot occur for other interaction than dipole-dipole

                // Angular dependent dipole-dipole interaction
                if (selectionRulesMultipole(r.state.first(), c.state.first(), 1) && selectionRulesMultipole(r.state.second(), c.state.second(), 1)) {
                    if (q1 == 0 && q2 == 0 && calculation_required[1]) {
                        scalar_t val = inverse_electric_constant * matrixelements1.getMultipole(r.state.first(), c.state.first(), 1) *
                                matrixelements2.getMultipole(r.state.second(), c.state.second(), 1);

                        this->addTriplet(interaction_angulardipole_triplets[1], r.idx, c.idx, val);

                    } else if (q1 != 0 && q2 != 0 && q1+q2 == 0 && (calculation_required[0] || calculation_required[2])) {
                        scalar_t val = inverse_electric_constant * matrixelements1.getMultipole(r.state.first(), c.state.first(), 1) *
                                matrixelements2.getMultipole(r.state.second(), c.state.second(), 1);

                        if (calculation_required[0]) this->addTriplet(interaction_angulardipole_triplets[0], r.idx, c.idx, val);
                        if (calculation_required[2]) this->addTriplet(interaction_angulardipole_triplets[2], r.idx, c.idx, -val);

                    } else if (std::abs(q1+q2) == 1 && calculation_required[3]) {
                        scalar_t val = inverse_electric_constant * matrixelements1.getMultipole(r.state.first(), c.state.first(), 1) *
                                matrixelements2.getMultipole(r.state.second(), c.state.second(), 1);

                        if (q1 == 1 || q2 == 1) this->addTriplet(interaction_angulardipole_triplets[3], r.idx, c.idx, -val);
                        else this->addTriplet(interaction_angulardipole_triplets[3], r.idx, c.idx, val);

                    } else if (std::abs(q1+q2) == 2 && calculation_required[2]) {
                        scalar_t val = inverse_electric_constant * matrixelements1.getMultipole(r.state.first(), c.state.first(), 1) *
                                matrixelements2.getMultipole(r.state.second(), c.state.second(), 1);

                        this->addTriplet(interaction_angulardipole_triplets[2], r.idx, c.idx, val);
                    }
                }

            } else {

                // Multipole interaction
                if (q1 == -q2) { // total momentum conserved
                    for (const auto &order : orange) {
                        double val = 0;
                        for (int kappa1 = 3; kappa1 <= order-2; ++kappa1) {
                            int kappa2 = order-1-kappa1;
                            if (selectionRulesMultipole(r.state.first(), c.state.first(), kappa1) && selectionRulesMultipole(r.state.second(), c.state.second(), kappa2)) {
                                double binomials = boost::math::binomial_coefficient<double>(kappa1+kappa2, kappa1+q1)*boost::math::binomial_coefficient<double>(kappa1+kappa2, kappa2-q2);
                                val += inverse_electric_constant * std::pow(-1,kappa2) * std::sqrt(binomials) * matrixelements1.getMultipole(r.state.first(), c.state.first(), kappa1)*
                                        matrixelements2.getMultipole(r.state.second(), c.state.second(), kappa2);
                            }
                        }

                        this->addTriplet(interaction_multipole_triplets[order], r.idx, c.idx, val);
                    }
                }

            }
        }
    }

    ////////////////////////////////////////////////////////////////////
    /// Build and transform the interaction to the used basis //////////
    ////////////////////////////////////////////////////////////////////

    for (size_t i = 0; i < calculation_required.size(); ++i) {
        if (!calculation_required[i]) continue;
        interaction_angulardipole[i].resize(states.size(),states.size());
        interaction_angulardipole[i].setFromTriplets(interaction_angulardipole_triplets[i].begin(), interaction_angulardipole_triplets[i].end());
        interaction_angulardipole_triplets[i].clear();

        interaction_angulardipole[i] = coefficients.adjoint()*interaction_angulardipole[i]*coefficients;
    }

    for (const auto &i : orange) {
        interaction_multipole[i].resize(states.size(),states.size());
        interaction_multipole[i].setFromTriplets(interaction_multipole_triplets[i].begin(), interaction_multipole_triplets[i].end());
        interaction_multipole_triplets[i].clear();

        interaction_multipole[i] = coefficients.adjoint()*interaction_multipole[i]*coefficients;
    }
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to construct Hamiltonian /////////
////////////////////////////////////////////////////////////////////

void SystemTwo::addInteraction() {
    if (distance == std::numeric_limits<double>::max()) return;

    // Build the total Hamiltonian
    double tolerance = 1e-24;

    if (angle != 0) { // setAngle and setOrder take care that a non-zero angle cannot occur for other interaction than dipole-dipole

        double powerlaw = 1./std::pow(distance, 3);
        for (size_t i = 0; i < 4; ++i) {
            if (std::abs(angle_terms[i]) > tolerance) hamiltonianmatrix += interaction_angulardipole[i]*angle_terms[i]*powerlaw;
        }

    } else {

        for (unsigned int order = 3; order <= ordermax; ++order) {
            double powerlaw = 1./std::pow(distance, order);
            hamiltonianmatrix += interaction_angulardipole[order]*powerlaw;
        }

    }
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to transform the interaction /////
////////////////////////////////////////////////////////////////////

void SystemTwo::transformInteraction(const eigen_sparse_t &transformator)  {
    for (auto &entry : interaction_angulardipole) entry.second = transformator.adjoint()*entry.second*transformator;
    for (auto &entry : interaction_multipole) entry.second = transformator.adjoint()*entry.second*transformator;
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to delete the interaction ////////
////////////////////////////////////////////////////////////////////

void SystemTwo::deleteInteraction()  {
    interaction_angulardipole.clear();
    interaction_multipole.clear();
}

////////////////////////////////////////////////////////////////////
/// Utility methods ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

void SystemTwo::addCoefficient(const size_t &row_1, const size_t &row_2, const size_t &col_new, const scalar_t &value_new, std::vector<eigen_triplet_t> &coefficients_triplets, std::vector<double> &sqnorm_list) {
    StateTwo state = StateTwo(system1.getStates()[row_1], system2.getStates()[row_2]);
    auto state_iter = states.get<1>().find(state);

    size_t row_new;
    if (state_iter != states.get<1>().end()) {
        row_new = state_iter->idx;
    } else {
        row_new = states.size();
        states.push_back(enumerated_state<StateTwo>(row_new, state));
    }

    coefficients_triplets.push_back(eigen_triplet_t(row_new, col_new, value_new));
    sqnorm_list[row_new] += std::pow(std::abs(value_new), 2);
}

void SystemTwo::addTriplet(std::vector<eigen_triplet_t> &triplets, const size_t r_idx, const size_t c_idx, const scalar_t val) {
    triplets.push_back(eigen_triplet_t(r_idx, c_idx, val));
    if (r_idx != c_idx) triplets.push_back(eigen_triplet_t(c_idx, r_idx, this->conjugate(val))); // triangular matrix is not sufficient because of basis change
}
