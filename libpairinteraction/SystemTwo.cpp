/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke, Johannes Block. All rights reserved.
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

#include "SystemTwo.h"
#include "GreenTensor.h"
#include "dtypes.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_set>
#include <vector>

SystemTwo::SystemTwo(const SystemOne &b1, const SystemOne &b2, MatrixElementCache &cache)
    : SystemBase(cache), species({{b1.getElement(), b2.getElement()}}), system1(b1), system2(b2),
      distance(std::numeric_limits<double>::max()), angle(0), ordermax(3), sym_permutation(NA),
      sym_inversion(NA), sym_reflection(NA), sym_rotation({ARB}) {}

SystemTwo::SystemTwo(const SystemOne &b1, const SystemOne &b2, MatrixElementCache &cache,
                     bool memory_saving)
    : SystemBase(cache, memory_saving), species({{b1.getElement(), b2.getElement()}}), system1(b1),
      system2(b2), distance(std::numeric_limits<double>::max()), angle(0), ordermax(3),
      sym_permutation(NA), sym_inversion(NA), sym_reflection(NA), sym_rotation({ARB}) {}

std::vector<StateOne>
SystemTwo::getStatesFirst() { // TODO @hmenke typemap for "state_set<StateOne>"
    this->buildBasis();
    std::unordered_set<StateOne> states_one;
    for (const auto &entry : states) {
        states_one.insert(StateOne(entry.state.getFirstState()));
    }
    return std::vector<StateOne>(states_one.begin(), states_one.end());
}

std::vector<StateOne>
SystemTwo::getStatesSecond() { // TODO @hmenke typemap for "state_set<StateOne>"
    this->buildBasis();
    std::unordered_set<StateOne> states_one;
    for (const auto &entry : states) {
        states_one.insert(StateOne(entry.state.getSecondState()));
    }
    return std::vector<StateOne>(states_one.begin(), states_one.end());
}

const std::array<std::string, 2> &SystemTwo::getElement() { return species; }

void SystemTwo::setDistance(double d) {
    this->onParameterChange();
    distance = d;
}

void SystemTwo::setDistanceX(double xAB) {
    this->onParameterChange();
    x = xAB;
}

void SystemTwo::setDistanceZA(double za) {
    this->onParameterChange();
    zA = za;
}
void SystemTwo::setDistanceZB(double zb) {
    this->onParameterChange();
    zB = zb;
}
void SystemTwo::setGTbool(bool GTboolean) {
    this->onParameterChange();
    GTbool = GTboolean;
}


void SystemTwo::setAngle(double a) {
    if (a != 0 && ordermax > 3) {
        throw std::runtime_error("A non-zero interaction angle can be directly used only for "
                                 "dipole-dipole interaction.");
    }

    this->onParameterChange();
    angle = a;

    angle_terms[0] = -1.;
    angle_terms[1] = 1. - 3. * std::pow(std::cos(angle), 2);
    angle_terms[2] = -1.5 * std::pow(std::sin(angle), 2);
    angle_terms[3] = -3. / std::sqrt(2) * std::sin(angle) * std::cos(angle);
}

void SystemTwo::setOrder(double o) {
    if (angle != 0 && o > 3) {
        throw std::runtime_error("A non-zero interaction angle can be directly used only for "
                                 "dipole-dipole interaction.");
    }

    this->onParameterChange();
    ordermax = o;
}

void SystemTwo::setConservedParityUnderPermutation(parity_t parity) {
    this->onSymmetryChange();
    sym_permutation = parity;
}

void SystemTwo::setConservedParityUnderInversion(parity_t parity) {
    this->onSymmetryChange();
    sym_inversion = parity;
}

void SystemTwo::setConservedParityUnderReflection(parity_t parity) {
    this->onSymmetryChange();
    sym_reflection = parity;
    if (!this->isRefelectionAndRotationCompatible()) {
        throw std::runtime_error("The conserved parity under reflection is not compatible to the "
                                 "previously specified conserved momenta.");
    }
    // if (sym_reflection != NA) std::cerr << "Warning: The one-atom states must already be
    // reflection symmetric in order to build reflection symmetric two-atom states." << std::endl; //
    // TODO make it work with one-atom states that are not pre-symmetrized
}

void SystemTwo::setConservedMomentaUnderRotation(const std::set<int> &momenta) {
    if (momenta.count(ARB) != 0 && momenta.size() > 1) {
        throw std::runtime_error(
            "If ARB (=arbitrary momentum) is specified, momenta must not be passed explicitely.");
    }
    this->onSymmetryChange();
    sym_rotation = momenta;
    if (!this->isRefelectionAndRotationCompatible()) {
        throw std::runtime_error("The conserved momenta are not compatible to the previously "
                                 "specified conserved parity under reflection.");
    }
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to initialize Basis //////////////
////////////////////////////////////////////////////////////////////

void SystemTwo::initializeBasis() {
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

    if (sym_permutation != NA || sym_permutation != NA) {
        // TODO check system1 == system2
    }

    ////////////////////////////////////////////////////////////////////
    /// Combine one atom states ////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // TODO consider further symmetries and check whether they are applicable

    std::vector<eigen_triplet_t> hamiltonianmatrix_triplets;
    hamiltonianmatrix_triplets.reserve(system1.getNumVectors() * system2.getNumVectors());
    states.reserve(system1.getNumStates() * system2.getNumStates());
    std::vector<eigen_triplet_t> coefficients_triplets; // TODO reserve
    std::vector<double> sqnorm_list(system1.getNumStates() * system2.getNumStates(), 0);

    size_t col_new = 0;
    for (size_t col_1 = 0; col_1 < system1.getNumVectors(); ++col_1) {
        for (size_t col_2 = 0; col_2 < system2.getNumVectors(); ++col_2) {

            // In case of inversion symmetry: skip half of the basis vector pairs
            if ((sym_inversion == EVEN && col_1 <= col_2) || // gerade
                (sym_inversion == ODD && col_1 < col_2)) {   // ungerade
                continue;
            }

            // In case of permutation symmetry: skip half of the basis vector pairs
            if ((sym_permutation == EVEN && col_1 <= col_2) || // sym
                (sym_permutation == ODD && col_1 < col_2)) {   // asym
                continue;
            }

            // Continue if the pair statet energy is not valid
            double energy = this->real(system1.getHamiltonianmatrix().coeff(col_1, col_1) +
                                       system2.getHamiltonianmatrix().coeff(col_2, col_2));
            if (!checkIsEnergyValid(energy)) {
                continue;
            }

            // Store the pair state energy
            hamiltonianmatrix_triplets.emplace_back(col_new, col_new, energy);

            // Build the basis vector that corresponds to the stored pair state energy
            for (eigen_iterator_t triple_1(system1.getCoefficients(), col_1); triple_1;
                 ++triple_1) {
                size_t row_1 = triple_1.row();
                StateOne state_1 = system1.getStates()[row_1]; // TODO cache states before

                for (eigen_iterator_t triple_2(system2.getCoefficients(), col_2); triple_2;
                     ++triple_2) {
                    size_t row_2 = triple_2.row();
                    StateOne state_2 = system2.getStates()[row_2]; // TODO cache states before

                    scalar_t value_new = triple_1.value() * triple_2.value();

                    int M = state_1.m + state_2.m;
                    int parityL = std::pow(-1, state_1.l + state_2.l);
                    int parityJ = std::pow(-1, state_1.j + state_2.j);
                    int parityM = std::pow(-1, M);

                    // Consider rotation symmetry
                    if (sym_rotation.count(ARB) == 0 && sym_rotation.count(M) == 0) {
                        continue;
                    }

                    // Combine symmetries
                    bool skip_reflection = false;
                    if (col_1 != col_2) {
                        // In case of inversion and permutation symmetry: the inversion symmetric
                        // state is already permutation symmetric
                        if (sym_inversion != NA && sym_permutation != NA) {
                            if (((sym_inversion == EVEN) ? -parityL : parityL) !=
                                ((sym_permutation == EVEN) ? -1 : 1)) {
                                continue; // the parity under inversion and permutation is different
                            }
                        }

                        // In case of inversion or permutation and reflection symmetry: the
                        // inversion or permutation symmetric state is already reflection symmetric
                        if ((sym_inversion != NA || sym_permutation != NA) &&
                            sym_reflection != NA &&
                            StateTwo({{state_1.species, state_2.species}}, {{state_1.n, state_2.n}},
                                     {{state_1.l, state_2.l}}, {{state_1.j, state_2.j}},
                                     {{-state_1.m, -state_2.m}}) == StateTwo(state_2, state_1)) {
                            if (sym_inversion != NA) {
                                if (((sym_inversion == EVEN) ? -parityL : parityL) !=
                                    ((sym_reflection == EVEN) ? parityL * parityJ * parityM
                                                              : -parityL * parityJ * parityM)) {
                                    continue; // the parity under inversion and reflection is
                                              // different
                                }
                                skip_reflection =
                                    true; // the parity under inversion and reflection is the same

                            } else if (sym_permutation != NA) {
                                if (((sym_permutation == EVEN) ? -1 : 1) !=
                                    ((sym_reflection == EVEN) ? parityL * parityJ * parityM
                                                              : -parityL * parityJ * parityM)) {
                                    continue; // the parity under permutation and reflection is
                                              // different
                                }
                                skip_reflection =
                                    true; // the parity under permutation and reflection is the same
                            }
                        }
                    }

                    // Adapt the normalization if required by symmetries
                    if (col_1 != col_2) {
                        if (sym_inversion != NA || sym_permutation != NA) {
                            value_new /= std::sqrt(2);
                        }
                    }
                    if (sym_reflection != NA && !skip_reflection) {
                        value_new /=
                            std::sqrt(2) *
                            std::sqrt(
                                2); // the second factor std::sqrt(2) is because of double counting
                    }

                    // Add an entry to the current basis vector
                    this->addCoefficient(StateTwo(state_1, state_2), col_new, value_new,
                                         coefficients_triplets, sqnorm_list);

                    // Add further entries to the current basis vector if required by symmetries
                    if (col_1 != col_2) {
                        if (sym_inversion != NA) {
                            scalar_t v = value_new;
                            v *= (sym_inversion == EVEN) ? -parityL : parityL;
                            this->addCoefficient(StateTwo(state_2, state_1), col_new, v,
                                                 coefficients_triplets, sqnorm_list);
                        } else if (sym_permutation != NA) {
                            scalar_t v = value_new;
                            v *= (sym_permutation == EVEN) ? -1 : 1;
                            this->addCoefficient(StateTwo(state_2, state_1), col_new, v,
                                                 coefficients_triplets, sqnorm_list);
                        }
                    }

                    if (sym_reflection != NA && !skip_reflection) {
                        scalar_t v = value_new;
                        v *= (sym_reflection == EVEN) ? parityL * parityJ * parityM
                                                      : -parityL * parityJ * parityM;
                        this->addCoefficient(
                            StateTwo({{state_1.species, state_2.species}}, {{state_1.n, state_2.n}},
                                     {{state_1.l, state_2.l}}, {{state_1.j, state_2.j}},
                                     {{-state_1.m, -state_2.m}}),
                            col_new, v, coefficients_triplets, sqnorm_list);

                        if (col_1 != col_2) {
                            if (sym_inversion != NA) {
                                scalar_t v = value_new;
                                v *= (sym_reflection == EVEN) ? parityL * parityJ * parityM
                                                              : -parityL * parityJ * parityM;
                                v *= (sym_inversion == EVEN) ? -parityL : parityL;
                                this->addCoefficient(
                                    StateTwo({{state_2.species, state_1.species}},
                                             {{state_2.n, state_1.n}}, {{state_2.l, state_1.l}},
                                             {{state_2.j, state_1.j}}, {{-state_2.m, -state_1.m}}),
                                    col_new, v, coefficients_triplets, sqnorm_list);
                            } else if (sym_permutation != NA) {
                                scalar_t v = value_new;
                                v *= (sym_reflection == EVEN) ? parityL * parityJ * parityM
                                                              : -parityL * parityJ * parityM;
                                v *= (sym_permutation == EVEN) ? -1 : 1;
                                this->addCoefficient(
                                    StateTwo({{state_2.species, state_1.species}},
                                             {{state_2.n, state_1.n}}, {{state_2.l, state_1.l}},
                                             {{state_2.j, state_1.j}}, {{-state_2.m, -state_1.m}}),
                                    col_new, v, coefficients_triplets, sqnorm_list);
                            }
                        }
                    }
                }
            }

            ++col_new;
        }
    }

    // Delete unecessary storage
    // system1 = SystemOne(species[0], cache); // TODO
    // system2 = SystemOne(species[1], cache); // TODO

    // Build data
    states.shrink_to_fit();

    coefficients.resize(states.size(), col_new);
    coefficients.setFromTriplets(coefficients_triplets.begin(), coefficients_triplets.end());
    coefficients_triplets.clear();

    hamiltonianmatrix.resize(col_new, col_new);
    hamiltonianmatrix.setFromTriplets(hamiltonianmatrix_triplets.begin(),
                                      hamiltonianmatrix_triplets.end());
    hamiltonianmatrix_triplets.clear();

    ////////////////////////////////////////////////////////////////////
    /// Remove vectors with too small norm /////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // Build transformator and remove vectors (if the squared norm is too small)
    std::vector<eigen_triplet_t> triplets_transformator;
    triplets_transformator.reserve(coefficients.cols());

    size_t idx_new = 0;
    for (int idx = 0; idx < coefficients.cols(); ++idx) { // idx = col = num basis vector
        double_t sqnorm = 0;

        // Calculate the square norm of the columns of the coefficient matrix
        for (eigen_iterator_t triple(coefficients, idx); triple; ++triple) {
            sqnorm += std::pow(std::abs(triple.value()), 2);
        }

        if (sqnorm > 0.05) {
            triplets_transformator.emplace_back(idx, idx_new++, 1);
        }
    }

    this->applyRightsideTransformator(triplets_transformator);

    ////////////////////////////////////////////////////////////////////
    /// Remove states that barely occur ////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // Build transformator and remove states if the squared norm is to small
    removeRestrictedStates([=](const enumerated_state<StateTwo> &entry) -> bool {
        return sqnorm_list[entry.idx] > 0.05;
    });
    // TODO CHECK Here for LeRoy-radius of largest quantum numbers
    /*////////////////////////////////////////////////////////////////////
    /// Sort states ////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // TODO put this into an extra method

    // Initialize original index locations
    std::vector<size_t> idx_sorted(states.size());
    iota(idx_sorted.begin(), idx_sorted.end(), 0);

    // Sort indexes based on comparing values in v
    std::sort(idx_sorted.begin(), idx_sorted.end(), [=](size_t i1, size_t i2) {return
    this->states[i1].state < this->states[i2].state;});

    // Make use of the sorted indexes in order to sort the states and transform the coefficients
    accordingly states_set<StateTwo> states_new; states_new.reserve(states.size());
    triplets_transformator.clear();
    triplets_transformator.reserve(states.size());

    idx_new = 0;
    for (size_t idx : idx_sorted) {
        states_new.push_back(enumerated_state<StateTwo>(idx_new, states[idx].state));
        triplets_transformator.push_back(eigen_triplet_t(idx_new++,idx,1));
    }

    states_new.shrink_to_fit();
    states = states_new;
    this->applyLeftsideTransformator(triplets_transformator);*/
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to calculate the interaction /////
////////////////////////////////////////////////////////////////////

void SystemTwo::initializeInteraction() {
    if (distance == std::numeric_limits<double>::max()) {
        return;
    }

    ////////////////////////////////////////////////////////////////////
    /// Prepare the calculation of the interaction /////////////////////
    ////////////////////////////////////////////////////////////////////

    // Check if something to do
    double tolerance = 1e-24;

    std::vector<bool> calculation_required(4, false);
    std::vector<int> orange;

    if (angle != 0) { // setAngle and setOrder take care that a non-zero angle cannot occur for
                      // other interaction than dipole-dipole

        for (size_t i = 0; i < 4; ++i) {
            if (std::abs(angle_terms[i]) > tolerance &&
                interaction_angulardipole.find(i) == interaction_angulardipole.end()) {
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
    auto states1 = this->getStatesFirst();
    auto states2 = this->getStatesSecond();

    for (unsigned int kappa = 1; kappa <= ordermax - 2; ++kappa) {
        cache.precalculateMultipole(states1, kappa);
        cache.precalculateMultipole(states2, kappa); // TODO check whether system1 == system2
    }
    ////////////////////////////////////////////////////////////////////
    /// Initialize dipolemoments as vectors first //////////////////////
    ////////////////////////////////////////////////////////////////////
    

    ////////////////////////////////////////////////////////////////////
    /// Generate the interaction in the canonical basis ////////////////
    ////////////////////////////////////////////////////////////////////

    std::unordered_map<int, std::vector<eigen_triplet_t>>
        interaction_angulardipole_triplets; // TODO reserve
    std::unordered_map<int, std::vector<eigen_triplet_t>>
        interaction_multipole_triplets; // TODO reserve
        
        
    if(GTbool){
        //matrix definitions in SystemTwo.h

      double dipolemoment1,dipolemoment2;
      std::complex<double> imagunit = std::complex<double> (0.,1.);
      std::complex<double> vec1[3];
      std::complex<double> vec2[3];
      std::complex<double> value;
      for (const auto &c : states) {
	for (const auto &r : states) {
	  
          dipolemoment1 = coulombs_constant*cache.getElectricDipole(r.state.first(), c.state.first());
          dipolemoment2 = coulombs_constant*cache.getElectricDipole(r.state.second(), c.state.second());     
            if(r.state.first().m - c.state.first().m == -1.){
                vec1[0] = (std::sqrt(1./2.)*(-dipolemoment1));
                vec1[1] = (std::sqrt(1./2.)*imagunit*dipolemoment1);
                vec1[2] = 0.;
            }
            if(r.state.second().m - c.state.second().m == -1.){
                vec2[0] = (std::sqrt(1./2.)*(-dipolemoment2));
                vec2[1] = (std::sqrt(1./2.)*imagunit*dipolemoment2);
                vec2[2] = 0.;
            }
            if(r.state.first().m - c.state.first().m == -1.){
                vec1[0] = (std::sqrt(1./2.)*(-dipolemoment1));
                vec1[1] = (std::sqrt(1./2.)*imagunit*dipolemoment1);
                vec1[2] = (0.);
            }
            if(r.state.second().m - c.state.second().m== -1.){
                vec2[0] = (std::sqrt(1./2.)*(-dipolemoment2));
                vec2[1] = (std::sqrt(1./2.)*imagunit*dipolemoment2);
                vec2[2] = (0.);
            }
            if(r.state.first().m - c.state.first().m == 0.){
                vec1[0] = 0.;
                vec1[1] = 0.;
                vec1[2] = std::sqrt(1.)*dipolemoment1;
            }
            if(r.state.second().m - c.state.second().m == 0.){
                vec2[0] = 0.;
                vec2[1] = 0.;
                vec2[2] = std::sqrt(1.)*dipolemoment2;
            }
            if(vec1[0]*vec2[0]!=0.){
                value = vec1[0]*vec2[0];
//                 value = vec1[0]*vec2[0]*GT.tensor(0,0); // => Dann für jeden GT alles neu ausrechenn. Das ist zu langwierig. Dieser Schritt wird in andere Funktion ausgelagert.
                this->addTripletC(xxGTmatrix,r.idx,c.idx,value);
            }
            if(vec1[1]*vec2[1]!=0.){
                value = vec1[1]*vec2[1];
                this->addTripletC(yyGTmatrix,r.idx,c.idx,value);
            }
            if(vec1[2]*vec2[2]!=0.){
                value = vec1[2]*vec2[2];
                this->addTripletC(zzGTmatrix,r.idx,c.idx,value);
            }
            if(vec1[0]*vec2[2]!=0.){
                value = vec1[0]*vec2[2];
                this->addTripletC(xzGTmatrix,r.idx,c.idx,value);
            }
            if(vec1[2]*vec2[0]!=0.){
                value = vec1[2]*vec2[0];
                this->addTripletC(zxGTmatrix,r.idx,c.idx,value);
            }
        }
      }
		
	    
    }
    
    for (const auto &c : states) { // TODO parallelization
        if (c.state.species.empty()) {
            continue; // TODO artifical states TODO [dummystates]
        }

        // Loop over row entries
        for (const auto &r : states) {
            if (r.state.species.empty()) {
                continue; // TODO artifical states TODO [dummystates]
            }
            if (r.idx < c.idx) {
                continue;
            }

            int q1 = r.state.first().m - c.state.first().m;
            int q2 = r.state.second().m - c.state.second().m;

            if (angle != 0) { // setAngle and setOrder take care that a non-zero angle cannot occur
                              // for other interaction than dipole-dipole

                // Angular dependent dipole-dipole interaction
                if (selectionRulesMultipoleNew(r.state.first(), c.state.first(), 1) &&
                    selectionRulesMultipoleNew(r.state.second(), c.state.second(), 1)) {
                    if (q1 == 0 && q2 == 0 && calculation_required[1]) {
                        scalar_t val = coulombs_constant *
                            cache.getElectricDipole(r.state.first(), c.state.first()) *
                            cache.getElectricDipole(r.state.second(), c.state.second());

                        this->addTriplet(interaction_angulardipole_triplets[1], r.idx, c.idx, val);

                    } else if (q1 != 0 && q2 != 0 && q1 + q2 == 0 &&
                               (calculation_required[0] || calculation_required[2])) {
                        scalar_t val = coulombs_constant *
                            cache.getElectricDipole(r.state.first(), c.state.first()) *
                            cache.getElectricDipole(r.state.second(), c.state.second());

                        if (calculation_required[0]) {
                            this->addTriplet(interaction_angulardipole_triplets[0], r.idx, c.idx,
                                             val);
                        }
                        if (calculation_required[2]) {
                            this->addTriplet(interaction_angulardipole_triplets[2], r.idx, c.idx,
                                             -val);
                        }

                    } else if (std::abs(q1 + q2) == 1 && calculation_required[3]) {
                        scalar_t val = coulombs_constant *
                            cache.getElectricDipole(r.state.first(), c.state.first()) *
                            cache.getElectricDipole(r.state.second(), c.state.second());

                        if (q1 == 1 || q2 == 1) {
                            this->addTriplet(interaction_angulardipole_triplets[3], r.idx, c.idx,
                                             -val);
                        } else {
                            this->addTriplet(interaction_angulardipole_triplets[3], r.idx, c.idx,
                                             val);
                        }

                    } else if (std::abs(q1 + q2) == 2 && calculation_required[2]) {
                        scalar_t val = coulombs_constant *
                            cache.getElectricDipole(r.state.first(), c.state.first()) *
                            cache.getElectricDipole(r.state.second(), c.state.second());

                        this->addTriplet(interaction_angulardipole_triplets[2], r.idx, c.idx, val);
                    }
                }

            } 
            
            else {
                // Multipole interaction
                if (q1 + q2 == 0) { // total momentum conserved
                    for (const auto &order : orange) {
                        double val = 0;
                        for (int kappa1 = 1; kappa1 <= order - 2; ++kappa1) {
                            int kappa2 = order - 1 - kappa1;
                            if (selectionRulesMultipoleNew(r.state.first(), c.state.first(),
                                                           kappa1) &&
                                selectionRulesMultipoleNew(r.state.second(), c.state.second(),
                                                           kappa2)) {
                                double binomials = boost::math::binomial_coefficient<double>(
                                                       kappa1 + kappa2, kappa1 + q1) *
                                    boost::math::binomial_coefficient<double>(kappa1 + kappa2,
                                                                              kappa2 - q2);
                                val += coulombs_constant * std::pow(-1, kappa2) *
                                    std::sqrt(binomials) *
                                    cache.getElectricMultipole(r.state.first(), c.state.first(),
                                                               kappa1) *
                                    cache.getElectricMultipole(r.state.second(), c.state.second(),
                                                               kappa2);
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
        if (!calculation_required[i]) {
            continue;
        }
        interaction_angulardipole[i].resize(states.size(), states.size());
        interaction_angulardipole[i].setFromTriplets(interaction_angulardipole_triplets[i].begin(),
                                                     interaction_angulardipole_triplets[i].end());
        interaction_angulardipole_triplets[i].clear();

        interaction_angulardipole[i] =
            coefficients.adjoint() * interaction_angulardipole[i] * coefficients;
    }

    for (const auto &i : orange) {
        interaction_multipole[i].resize(states.size(), states.size());
        interaction_multipole[i].setFromTriplets(interaction_multipole_triplets[i].begin(),
                                                 interaction_multipole_triplets[i].end());
        interaction_multipole_triplets[i].clear();

        interaction_multipole[i] = coefficients.adjoint() * interaction_multipole[i] * coefficients;
    }
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to construct Hamiltonian /////////
////////////////////////////////////////////////////////////////////

void SystemTwo::addInteraction() {
    if (distance == std::numeric_limits<double>::max()) {
        return;
    }

    // Build the total Hamiltonian
    double tolerance = 1e-24;
    
    if(GTbool){
      //Ganz umständlich?:
      GreenTensor GT(x,zA,zB);
      Eigen::SparseMatrix<std::complex<double> > dummy;
      Eigen::SparseMatrix<std::complex<double> > GThamiltonian;
      dummy.setFromTriplets(xxGTmatrix.begin(),xxGTmatrix.end());
      dummy = dummy*GT.tensor(0,0);
      GThamiltonian = dummy;
      dummy.setFromTriplets(yyGTmatrix.begin(),yyGTmatrix.end());
      dummy = dummy*GT.tensor(1,1);
      GThamiltonian += dummy;
      dummy.setFromTriplets(zzGTmatrix.begin(),zzGTmatrix.end());
      dummy = dummy*GT.tensor(2,2);
      GThamiltonian += dummy;
      dummy.setFromTriplets(xzGTmatrix.begin(),xzGTmatrix.end());
      dummy = dummy*GT.tensor(0,2);
      GThamiltonian += dummy;
      dummy.setFromTriplets(zxGTmatrix.begin(),zxGTmatrix.end());
      dummy = dummy*GT.tensor(2,0);
      GThamiltonian += dummy;
    }



    else if (angle != 0) { // setAngle and setOrder take care that a non-zero angle cannot occur for
                      // other interaction than dipole-dipole

        double powerlaw = 1. / std::pow(distance, 3);
        for (size_t i = 0; i < 4; ++i) {
            if (std::abs(angle_terms[i]) > tolerance) {
                hamiltonianmatrix += interaction_angulardipole[i] * angle_terms[i] * powerlaw;
            }
        }

    } else {

        for (unsigned int order = 3; order <= ordermax; ++order) {
            double powerlaw = 1. / std::pow(distance, order);
            hamiltonianmatrix += interaction_multipole[order] * powerlaw;
        }
    }
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to construct Hamiltonian using GreenTensor /////////
////////////////////////////////////////////////////////////////////



void SystemTwo::DipoleVector(){
    //just empty, everything moved to other place. Leave it as place holder for now.
    
}


////////////////////////////////////////////////////////////////////
/// Method that allows base class to transform the interaction /////
////////////////////////////////////////////////////////////////////

void SystemTwo::transformInteraction(const eigen_sparse_t &transformator) {
    for (auto &entry : interaction_angulardipole) {
        entry.second = transformator.adjoint() * entry.second * transformator;
    }
    for (auto &entry : interaction_multipole) {
        entry.second = transformator.adjoint() * entry.second * transformator; // NOLINT
    }
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to delete the interaction ////////
////////////////////////////////////////////////////////////////////

void SystemTwo::deleteInteraction() {
    interaction_angulardipole.clear();
    interaction_multipole.clear();
}

////////////////////////////////////////////////////////////////////
/// Methods that allows base class to rotate states ////////////////
////////////////////////////////////////////////////////////////////

eigen_sparse_t SystemTwo::rotateStates(const std::vector<size_t> &states_indices, double alpha,
                                       double beta, double gamma) {
    // Initialize Wigner D matrix
    WignerD wigner;

    // Rotate state
    std::vector<eigen_triplet_t> states_rotated_triplets;
    states_rotated_triplets.reserve(
        std::min(static_cast<size_t>(10 * 10), states.size()) *
        states_indices.size()); // TODO std::min( (2*jmax+1)*(2*jmax+1), states.size() ) *
                                // states_indices.size()

    size_t current = 0;
    for (auto const &idx : states_indices) {
        this->addRotated(states[idx].state, current++, states_rotated_triplets, wigner, alpha, beta,
                         gamma);
    }

    eigen_sparse_t states_rotated(states.size(), states_indices.size());
    states_rotated.setFromTriplets(states_rotated_triplets.begin(), states_rotated_triplets.end());
    states_rotated_triplets.clear();

    return states_rotated;
}

eigen_sparse_t SystemTwo::buildStaterotator(double alpha, double beta, double gamma) {
    // Initialize Wigner D matrix
    WignerD wigner;

    // Build rotator
    std::vector<eigen_triplet_t> rotator_triplets;
    rotator_triplets.reserve(
        std::min(static_cast<size_t>(10 * 10), states.size()) *
        states.size()); // TODO std::min( (2*jmax+1)*(2*jmax+1), states.size() ) * states.size()

    for (auto const &entry : states) {
        this->addRotated(entry.state, entry.idx, rotator_triplets, wigner, alpha, beta, gamma);
    }

    eigen_sparse_t rotator(states.size(), states.size());
    rotator.setFromTriplets(rotator_triplets.begin(), rotator_triplets.end()); // NOLINT
    rotator_triplets.clear();

    return rotator;
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to combine systems ///////////////
////////////////////////////////////////////////////////////////////

void SystemTwo::incorporate(SystemBase<StateTwo> &system) {
    // Combine parameters
    if (species[0] != dynamic_cast<SystemTwo &>(system).species[0]) {
        throw std::runtime_error(
            "The value of the variable 'element' must be the same for both systems.");
    }
    if (species[1] != dynamic_cast<SystemTwo &>(system).species[1]) {
        throw std::runtime_error(
            "The value of the variable 'element' must be the same for both systems.");
    }
    if (distance != dynamic_cast<SystemTwo &>(system).distance) {
        throw std::runtime_error(
            "The value of the variable 'distance' must be the same for both systems.");
    }
    if (angle != dynamic_cast<SystemTwo &>(system).angle) {
        throw std::runtime_error(
            "The value of the variable 'angle' must be the same for both systems."); // implies that
                                                                                     // angle_terms
                                                                                     // is the same,
                                                                                     // too
    }
    if (ordermax != dynamic_cast<SystemTwo &>(system).ordermax) {
        throw std::runtime_error(
            "The value of the variable 'ordermax' must be the same for both systems.");
    }

    // Combine symmetries
    unsigned int num_different_symmetries = 0;
    if (sym_permutation != dynamic_cast<SystemTwo &>(system).sym_permutation) {
        sym_permutation = NA;
        ++num_different_symmetries;
    }
    if (sym_inversion != dynamic_cast<SystemTwo &>(system).sym_inversion) {
        sym_inversion = NA;
        ++num_different_symmetries;
    }
    if (sym_reflection != dynamic_cast<SystemTwo &>(system).sym_reflection) {
        sym_reflection = NA;
        ++num_different_symmetries;
    }
    if (!std::equal(sym_rotation.begin(), sym_rotation.end(),
                    dynamic_cast<SystemTwo &>(system).sym_rotation.begin())) {
        if (sym_rotation.count(ARB) != 0 ||
            dynamic_cast<SystemTwo &>(system).sym_rotation.count(ARB) != 0) {
            sym_rotation = {ARB};
        } else {
            sym_rotation.insert(dynamic_cast<SystemTwo &>(system).sym_rotation.begin(),
                                dynamic_cast<SystemTwo &>(system).sym_rotation.end());
        }
        ++num_different_symmetries;
    }
    // if (num_different_symmetries > 1) std::cerr << "Warning: The systems differ in more than one
    // symmetry. For the combined system, the notion of symmetries might be meaningless." <<
    // std::endl; // TODO let the user enable/disable this warning

    // Clear cached interaction
    this->deleteInteraction();
}

////////////////////////////////////////////////////////////////////
/// Utility methods ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

void SystemTwo::addCoefficient(const StateTwo &state, const size_t &col_new,
                               const scalar_t &value_new,
                               std::vector<eigen_triplet_t> &coefficients_triplets,
                               std::vector<double> &sqnorm_list) {
    auto state_iter = states.get<1>().find(state);

    size_t row_new;
    if (state_iter != states.get<1>().end()) {
        row_new = state_iter->idx;
    } else {
        row_new = states.size();
        states.push_back(enumerated_state<StateTwo>(row_new, state));
    }

    coefficients_triplets.emplace_back(row_new, col_new, value_new);
    sqnorm_list[row_new] += std::pow(std::abs(value_new), 2);
}

void SystemTwo::addTriplet(std::vector<eigen_triplet_t> &triplets, const size_t r_idx,
                           const size_t c_idx, const scalar_t val) {
    triplets.emplace_back(r_idx, c_idx, val);
    if (r_idx != c_idx) {
        triplets.emplace_back(
            c_idx, r_idx,
            this->conjugate(val)); // triangular matrix is not sufficient because of basis change
    }
}

void SystemTwo::addTripletC(std::vector<eigen_triplet_complex_t> &triplets, const size_t r_idx,
                           const size_t c_idx, const std::complex<double> val) {
    triplets.emplace_back(r_idx, c_idx, val);
    if (r_idx != c_idx) {
        triplets.emplace_back(
            c_idx, r_idx,
            this->conjugate(val)); // triangular matrix is not sufficient because of basis change
    }
}

template <>
double SystemTwo::convert(const std::complex<double> &val) {
    return val.real();
}

bool SystemTwo::isRefelectionAndRotationCompatible() {
    if (sym_rotation.count(ARB) != 0 || sym_reflection == NA) {
        return true;
    }

    for (const auto &s : sym_rotation) {
        if (sym_rotation.count(-s) == 0) {
            return false;
        }
    }

    return true;
}
