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

#include "SystemTwo.h"
#include "dtypes.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_set>
#include <vector>

SystemTwo::SystemTwo(const SystemOne &b1, const SystemOne &b2, MatrixElementCache &cache)
    : SystemBase(cache), species({{b1.getSpecies(), b2.getSpecies()}}), system1(b1), system2(b2),
      minimal_le_roy_radius(std::numeric_limits<double>::max()),
      distance(std::numeric_limits<double>::max()), angle(0), ordermax(3), sym_permutation(NA),
      sym_inversion(NA), sym_reflection(NA), sym_rotation({ARB}) {}

SystemTwo::SystemTwo(const SystemOne &b1, const SystemOne &b2, MatrixElementCache &cache,
                     bool memory_saving)
    : SystemBase(cache, memory_saving), species({{b1.getSpecies(), b2.getSpecies()}}), system1(b1),
      system2(b2), minimal_le_roy_radius(std::numeric_limits<double>::max()),
      distance(std::numeric_limits<double>::max()), angle(0), ordermax(3), sym_permutation(NA),
      sym_inversion(NA), sym_reflection(NA), sym_rotation({ARB}) {}

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

const std::array<std::string, 2> &SystemTwo::getSpecies() { return species; }

void SystemTwo::setDistance(double d) {
    this->onParameterChange();
    distance = d;
    this->checkDistance(d);
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
    // reflection symmetric in order to build reflection symmetric two-atom states." << std::endl;

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
    /// Check whether the single atom states fit to the symmetries /////
    ////////////////////////////////////////////////////////////////////

    if (sym_permutation != NA || sym_permutation != NA) {
        // TODO check system1 == system2
    }

    // TODO consider further symmetries and check whether they are applicable

    ////////////////////////////////////////////////////////////////////
    /// Check which basis vectors contain artificial states ////////////
    ////////////////////////////////////////////////////////////////////

    // TODO check whether system1 == system2

    std::vector<bool> artificial1(system1.getNumBasisvectors(), false);
    for (size_t col = 0; col < system1.getNumBasisvectors(); ++col) {
        for (eigen_iterator_t triple(system1.getBasisvectors(), col); triple; ++triple) {
            if (system1.getStatesMultiIndex()[triple.row()].state.isArtificial()) {
                artificial1[triple.col()] = true;
            }
        }
    }

    std::vector<bool> artificial2(system2.getNumBasisvectors(), false);
    for (size_t col = 0; col < system2.getNumBasisvectors(); ++col) {
        for (eigen_iterator_t triple(system2.getBasisvectors(), col); triple; ++triple) {
            if (system2.getStatesMultiIndex()[triple.row()].state.isArtificial()) {
                artificial2[triple.col()] = true;
            }
        }
    }

    ////////////////////////////////////////////////////////////////////
    /// Build two atom states //////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    /// Combine one atom states ////////////////////////////////////////

    std::vector<eigen_triplet_t> hamiltonianmatrix_triplets;
    hamiltonianmatrix_triplets.reserve(system1.getNumBasisvectors() * system2.getNumBasisvectors());
    states.reserve(system1.getNumStates() * system2.getNumStates() + states_to_add.size());
    std::vector<eigen_triplet_t> basisvectors_triplets; // TODO reserve
    std::vector<double> sqnorm_list(
        system1.getNumStates() * system2.getNumStates() + states_to_add.size(), 0);

    int M = 0;
    int parityL = 0;
    int parityJ = 0;
    int parityM = 0;

    size_t col_new = 0;
    for (size_t col_1 = 0; col_1 < system1.getNumBasisvectors(); ++col_1) {
        for (size_t col_2 = 0; col_2 < system2.getNumBasisvectors(); ++col_2) {

            // In case of artificial states, some symmetries won't work
            auto sym_inversion_local = sym_inversion;
            auto sym_reflection_local = sym_reflection;
            auto sym_rotation_local = sym_rotation;
            if (artificial1[col_1] || artificial1[col_2]) {
                if (sym_inversion_local != NA || sym_reflection_local != NA ||
                    sym_rotation_local.count(ARB) == 0) {
                    std::cerr
                        << "WARNING: Only permutation symmetry can be applied to artificial states."
                        << std::endl;
                }
                sym_inversion_local = NA;
                sym_reflection_local = NA;
                sym_rotation_local = std::set<int>({ARB});
            }

            // In case of inversion or permutation symmetry: skip half of the basis vector pairs
            if ((sym_inversion_local == EVEN && col_1 <= col_2) || // gerade
                (sym_inversion_local == ODD && col_1 < col_2)) {   // ungerade
                continue;
            }
            if ((sym_permutation == EVEN && col_1 <= col_2) || // sym
                (sym_permutation == ODD && col_1 < col_2)) {   // asym
                continue;
            }

            // Continue if the pair statet energy is not valid
            double energy = this->real(system1.getHamiltonian().coeff(col_1, col_1) +
                                       system2.getHamiltonian().coeff(col_2, col_2));
            if (!checkIsEnergyValid(energy)) {
                continue;
            }

            // Store the pair state energy
            hamiltonianmatrix_triplets.emplace_back(col_new, col_new, energy);

            // Build the basis vector that corresponds to the stored pair state energy
            for (eigen_iterator_t triple_1(system1.getBasisvectors(), col_1); triple_1;
                 ++triple_1) {
                size_t row_1 = triple_1.row();
                StateOne state_1 = system1.getStatesMultiIndex()[row_1].state;

                for (eigen_iterator_t triple_2(system2.getBasisvectors(), col_2); triple_2;
                     ++triple_2) {
                    size_t row_2 = triple_2.row();
                    StateOne state_2 = system2.getStatesMultiIndex()[row_2].state;

                    scalar_t value_new = triple_1.value() * triple_2.value();

                    if (!artificial1[col_1] && !artificial1[col_2]) {
                        M = state_1.getM() + state_2.getM();
                        parityL = std::pow(-1, state_1.getL() + state_2.getL());
                        parityJ = std::pow(-1, state_1.getJ() + state_2.getJ());
                        parityM = std::pow(-1, M);
                    }

                    bool different = col_1 != col_2;

                    // Consider rotation symmetry
                    if (sym_rotation_local.count(ARB) == 0 && sym_rotation_local.count(M) == 0) {
                        continue;
                    }

                    // Combine symmetries
                    bool skip_reflection = false;
                    if (different) {
                        // In case of inversion and permutation symmetry: the inversion symmetric
                        // state is already permutation symmetric
                        if (sym_inversion_local != NA && sym_permutation != NA) {
                            if (((sym_inversion_local == EVEN) ? -parityL : parityL) !=
                                ((sym_permutation == EVEN) ? -1 : 1)) {
                                continue; // parity under inversion and permutation is different
                            }
                        }

                        // In case of inversion or permutation and reflection symmetry: the
                        // inversion or permutation symmetric state is already reflection symmetric
                        if (sym_inversion_local != NA && sym_reflection_local != NA &&
                            StateTwo(state_1.getReflected(), state_2.getReflected()) ==
                                StateTwo(state_2, state_1)) {
                            if (((sym_inversion_local == EVEN) ? -parityL : parityL) !=
                                ((sym_reflection_local == EVEN) ? parityL * parityJ * parityM
                                                                : -parityL * parityJ * parityM)) {
                                continue; // parity under inversion and reflection is different
                            }
                            skip_reflection =
                                true; // parity under inversion and reflection is the same

                        } else if (sym_permutation != NA && sym_reflection_local != NA &&
                                   StateTwo(state_1.getReflected(), state_2.getReflected()) ==
                                       StateTwo(state_2, state_1)) {
                            if (((sym_permutation == EVEN) ? -1 : 1) !=
                                ((sym_reflection_local == EVEN) ? parityL * parityJ * parityM
                                                                : -parityL * parityJ * parityM)) {
                                continue; // parity under permutation and reflection is different
                            }
                            skip_reflection =
                                true; // parity under permutation and reflection is the same
                        }
                    }

                    // Adapt the normalization if required by symmetries
                    if ((sym_inversion_local != NA || sym_permutation != NA) && different) {
                        value_new /= std::sqrt(2);
                    }
                    if (sym_reflection_local != NA && !skip_reflection) {
                        value_new /= std::sqrt(2) * std::sqrt(2);
                        // the second factor std::sqrt(2) is because of double counting
                    }

                    // Add an entry to the current basis vector
                    this->addBasisvectors(StateTwo(state_1, state_2), col_new, value_new,
                                          basisvectors_triplets, sqnorm_list);

                    // Add further entries to the current basis vector if required by symmetries
                    if (different) {
                        if (sym_inversion_local != NA) {
                            scalar_t v = value_new;
                            v *= (sym_inversion_local == EVEN) ? -parityL : parityL;
                            this->addBasisvectors(StateTwo(state_2, state_1), col_new, v,
                                                  basisvectors_triplets, sqnorm_list);
                        } else if (sym_permutation != NA) {
                            scalar_t v = value_new;
                            v *= (sym_permutation == EVEN) ? -1 : 1;
                            this->addBasisvectors(StateTwo(state_2, state_1), col_new, v,
                                                  basisvectors_triplets, sqnorm_list);
                        }
                    }
                    if (sym_reflection_local != NA && !skip_reflection) {
                        scalar_t v = value_new;
                        v *= (sym_reflection_local == EVEN) ? parityL * parityJ * parityM
                                                            : -parityL * parityJ * parityM;
                        this->addBasisvectors(
                            StateTwo(state_1.getReflected(), state_2.getReflected()), col_new, v,
                            basisvectors_triplets, sqnorm_list);

                        if (different) {
                            if (sym_inversion_local != NA) {
                                scalar_t v = value_new;
                                v *= (sym_reflection_local == EVEN) ? parityL * parityJ * parityM
                                                                    : -parityL * parityJ * parityM;
                                v *= (sym_inversion_local == EVEN) ? -parityL : parityL;
                                this->addBasisvectors(
                                    StateTwo(state_2.getReflected(), state_1.getReflected()),
                                    col_new, v, basisvectors_triplets, sqnorm_list);
                            } else if (sym_permutation != NA) {
                                scalar_t v = value_new;
                                v *= (sym_reflection_local == EVEN) ? parityL * parityJ * parityM
                                                                    : -parityL * parityJ * parityM;
                                v *= (sym_permutation == EVEN) ? -1 : 1;
                                this->addBasisvectors(
                                    StateTwo(state_2.getReflected(), state_1.getReflected()),
                                    col_new, v, basisvectors_triplets, sqnorm_list);
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

    /// Loop over user-defined states //////////////////////////////////

    // Check that the user-defined states are not already contained in the list of states
    for (const auto &state : states_to_add) {
        if (states.get<1>().find(state) != states.get<1>().end()) {
            throw std::runtime_error("The state " + state.str() +
                                     " is already contained in the list of states.");
        }
        for (int idx = 0; idx < 2; ++idx) {
            if (!state.isArtificial(idx) && state.getSpecies(idx) != species[idx]) {
                throw std::runtime_error("The state " + state.str() + " is of the wrong species.");
            }
        }
    }

    // Warn if reflection symmetry is selected
    if (!states_to_add.empty() && sym_reflection != NA) {
        std::cerr << "WARNING: Reflection symmetry cannot be handled for user-defined states."
                  << std::endl;
    }

    // Add user-defined states
    for (const auto &state : states_to_add) {
        bool different = state.getFirstState() != state.getSecondState();

        // Get energy of the state
        double energy = 0;
        for (int idx = 0; idx < 2; ++idx) {
            if (!state.isArtificial(idx)) {
                energy += state.getEnergy(idx);
            }
        }

        // In case of artificial states, some symmetries won't work
        auto sym_inversion_local = sym_inversion;
        auto sym_rotation_local = sym_rotation;
        if (state.isArtificial(0) || state.isArtificial(1)) {
            if (sym_inversion_local != NA || sym_rotation_local.count(ARB) == 0) {
                std::cerr
                    << "WARNING: Only permutation symmetry can be applied to artificial states."
                    << std::endl;
            }
            sym_inversion_local = NA;
            sym_rotation_local = std::set<int>({ARB});
        } else {
            M = state.getM(0) + state.getM(1);
            parityL = std::pow(-1, state.getL(0) + state.getL(1));
        }

        // Consider rotation symmetry
        if (sym_rotation_local.count(ARB) == 0 && sym_rotation_local.count(M) == 0) {
            continue;
        }

        // Combine symmetries (in case of inversion and permutation symmetry: the inversion
        // symmetric state is already permutation symmetric)
        if (sym_inversion_local != NA && sym_permutation != NA && different) {
            if (((sym_inversion_local == EVEN) ? -parityL : parityL) !=
                ((sym_permutation == EVEN) ? -1 : 1)) {
                continue; // parity under inversion and permutation is different
            }
        }

        // Check whether the symmetries can be realized with the states available
        if ((sym_inversion_local != NA || sym_permutation != NA) && different) {
            auto state_changed = StateTwo(state.getSecondState(), state.getFirstState());
            if (states_to_add.find(state_changed) == states_to_add.end()) {
                throw std::runtime_error("The state " + state_changed.str() +
                                         " required by symmetries cannot be found.");
            }
        }

        // In case of inversion or permutation symmetry: skip half of the states
        if ((sym_inversion_local == EVEN &&
             state.getFirstState() <= state.getSecondState()) || // gerade
            (sym_inversion_local == ODD &&
             state.getFirstState() < state.getSecondState())) { // ungerade
            continue;
        }
        if ((sym_permutation == EVEN && state.getFirstState() <= state.getSecondState()) || // sym
            (sym_permutation == ODD && state.getFirstState() < state.getSecondState())) {   // asym
            continue;
        }

        // Store the energy of the two atom state
        hamiltonianmatrix_triplets.emplace_back(col_new, col_new, energy);

        // Adapt the normalization if required by symmetries
        scalar_t value_new = 1;
        if ((sym_inversion_local != NA || sym_permutation != NA) && different) {
            value_new /= std::sqrt(2);
        }

        // Add an entry to the current basis vector
        this->addBasisvectors(state, col_new, value_new, basisvectors_triplets, sqnorm_list);

        // Add further entries to the current basis vector if required by symmetries
        if (different) {
            if (sym_inversion_local != NA) {
                scalar_t v = value_new;
                v *= (sym_inversion_local == EVEN) ? -parityL : parityL;
                this->addBasisvectors(StateTwo(state.getSecondState(), state.getFirstState()),
                                      col_new, v, basisvectors_triplets, sqnorm_list);
            } else if (sym_permutation != NA) {
                scalar_t v = value_new;
                v *= (sym_permutation == EVEN) ? -1 : 1;
                this->addBasisvectors(StateTwo(state.getSecondState(), state.getFirstState()),
                                      col_new, v, basisvectors_triplets, sqnorm_list);
            }
        }

        ++col_new;
    }

    /// Build data /////////////////////////////////////////////////////

    states.shrink_to_fit();

    basisvectors.resize(states.size(), col_new);
    basisvectors.setFromTriplets(basisvectors_triplets.begin(), basisvectors_triplets.end());
    basisvectors_triplets.clear();

    hamiltonianmatrix.resize(col_new, col_new);
    hamiltonianmatrix.setFromTriplets(hamiltonianmatrix_triplets.begin(),
                                      hamiltonianmatrix_triplets.end());
    hamiltonianmatrix_triplets.clear();

    ////////////////////////////////////////////////////////////////////
    /// Remove vectors with too small norm /////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // Build transformator and remove vectors (if the squared norm is too small)
    std::vector<eigen_triplet_t> triplets_transformator;
    triplets_transformator.reserve(basisvectors.cols());

    size_t idx_new = 0;
    for (int idx = 0; idx < basisvectors.cols(); ++idx) { // idx = col = num basis vector
        double_t sqnorm = 0;

        // Calculate the square norm of the columns of the coefficient matrix
        for (eigen_iterator_t triple(basisvectors, idx); triple; ++triple) {
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

    ////////////////////////////////////////////////////////////////////
    /// Get minimal Le Roy radius //////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // Estimate minimal Le Roy radius
    StateTwo crucial_state{{{"None", "None"}}};
    for (const auto &e : states) {
        if (e.state.isArtificial(0) || e.state.isArtificial(1)) {
            continue;
        }

        auto n = e.state.getNStar();
        auto l = e.state.getL();

        double le_roy_radius = 2 * au2um *
            (std::sqrt(0.5 * n[0] * n[0] * (5 * n[0] * n[0] + 1 - 3 * l[0] * (l[0] + 1))) +
             std::sqrt(0.5 * n[1] * n[1] * (5 * n[1] * n[1] + 1 - 3 * l[1] * (l[1] + 1))));

        if (le_roy_radius < minimal_le_roy_radius) {
            minimal_le_roy_radius = le_roy_radius;
            crucial_state = e.state;
        }
    }

    // Calculate minimal Le Roy radius precisely
    if (crucial_state.isArtificial(0) || crucial_state.isArtificial(1)) {
        minimal_le_roy_radius = 0;
    } else {
        minimal_le_roy_radius = 2 *
            (std::sqrt(
                 cache.getRadial(crucial_state.getFirstState(), crucial_state.getFirstState(), 2)) +
             std::sqrt(cache.getRadial(crucial_state.getSecondState(),
                                       crucial_state.getSecondState(), 2)));
    }

    // Check whether distances are larger than the minimal Le Roy radius
    this->checkDistance(distance);

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

    // Make use of the sorted indexes in order to sort the states and transform the basisvectors
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
    /// Generate the interaction in the canonical basis ////////////////
    ////////////////////////////////////////////////////////////////////

    std::unordered_map<int, std::vector<eigen_triplet_t>>
        interaction_angulardipole_triplets; // TODO reserve
    std::unordered_map<int, std::vector<eigen_triplet_t>>
        interaction_multipole_triplets; // TODO reserve

    // Loop over column entries
    for (const auto &c : states) { // TODO parallelization
        if (c.state.isArtificial(0) || c.state.isArtificial(1)) {
            continue;
        }

        // Loop over row entries
        for (const auto &r : states) {
            if (r.state.isArtificial(0) || r.state.isArtificial(1)) {
                continue;
            }

            if (r.idx < c.idx) {
                continue;
            }

            int q1 = r.state.getM(0) - c.state.getM(0);
            int q2 = r.state.getM(1) - c.state.getM(1);

            if (angle != 0) { // setAngle and setOrder take care that a non-zero angle cannot occur
                              // for other interaction than dipole-dipole

                // Angular dependent dipole-dipole interaction
                if (selectionRulesMultipoleNew(r.state.getFirstState(), c.state.getFirstState(),
                                               1) &&
                    selectionRulesMultipoleNew(r.state.getSecondState(), c.state.getSecondState(),
                                               1)) {
                    if (q1 == 0 && q2 == 0 && calculation_required[1]) {
                        scalar_t val = coulombs_constant *
                            cache.getElectricDipole(r.state.getFirstState(),
                                                    c.state.getFirstState()) *
                            cache.getElectricDipole(r.state.getSecondState(),
                                                    c.state.getSecondState());

                        this->addTriplet(interaction_angulardipole_triplets[1], r.idx, c.idx, val);

                    } else if (q1 != 0 && q2 != 0 && q1 + q2 == 0 &&
                               (calculation_required[0] || calculation_required[2])) {
                        scalar_t val = coulombs_constant *
                            cache.getElectricDipole(r.state.getFirstState(),
                                                    c.state.getFirstState()) *
                            cache.getElectricDipole(r.state.getSecondState(),
                                                    c.state.getSecondState());

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
                            cache.getElectricDipole(r.state.getFirstState(),
                                                    c.state.getFirstState()) *
                            cache.getElectricDipole(r.state.getSecondState(),
                                                    c.state.getSecondState());

                        if (q1 == 1 || q2 == 1) {
                            this->addTriplet(interaction_angulardipole_triplets[3], r.idx, c.idx,
                                             -val);
                        } else {
                            this->addTriplet(interaction_angulardipole_triplets[3], r.idx, c.idx,
                                             val);
                        }

                    } else if (std::abs(q1 + q2) == 2 && calculation_required[2]) {
                        scalar_t val = coulombs_constant *
                            cache.getElectricDipole(r.state.getFirstState(),
                                                    c.state.getFirstState()) *
                            cache.getElectricDipole(r.state.getSecondState(),
                                                    c.state.getSecondState());

                        this->addTriplet(interaction_angulardipole_triplets[2], r.idx, c.idx, val);
                    }
                }

            } else {

                // Multipole interaction
                if (q1 + q2 == 0) { // total momentum conserved
                    for (const auto &order : orange) {
                        double val = 0;
                        for (int kappa1 = 1; kappa1 <= order - 2; ++kappa1) {
                            int kappa2 = order - 1 - kappa1;
                            if (selectionRulesMultipoleNew(r.state.getFirstState(),
                                                           c.state.getFirstState(), kappa1) &&
                                selectionRulesMultipoleNew(r.state.getSecondState(),
                                                           c.state.getSecondState(), kappa2)) {
                                double binomials = boost::math::binomial_coefficient<double>(
                                                       kappa1 + kappa2, kappa1 + q1) *
                                    boost::math::binomial_coefficient<double>(kappa1 + kappa2,
                                                                              kappa2 - q2);
                                val += coulombs_constant * std::pow(-1, kappa2) *
                                    std::sqrt(binomials) *
                                    cache.getElectricMultipole(r.state.getFirstState(),
                                                               c.state.getFirstState(), kappa1) *
                                    cache.getElectricMultipole(r.state.getSecondState(),
                                                               c.state.getSecondState(), kappa2);
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
            basisvectors.adjoint() * interaction_angulardipole[i] * basisvectors;
    }

    for (const auto &i : orange) {
        interaction_multipole[i].resize(states.size(), states.size());
        interaction_multipole[i].setFromTriplets(interaction_multipole_triplets[i].begin(),
                                                 interaction_multipole_triplets[i].end());
        interaction_multipole_triplets[i].clear();

        interaction_multipole[i] = basisvectors.adjoint() * interaction_multipole[i] * basisvectors;
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

    if (angle != 0) { // setAngle and setOrder take care that a non-zero angle cannot occur for
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

void SystemTwo::checkDistance(double distance) {
    if (minimal_le_roy_radius != std::numeric_limits<double>::max() &&
        distance < minimal_le_roy_radius) {
        std::cerr << "WARNING: The distance " << distance
                  << " um is smaller than the Le Roy radius " << minimal_le_roy_radius << " um."
                  << std::endl;
    }
}

void SystemTwo::addBasisvectors(const StateTwo &state, const size_t &col_new,
                                const scalar_t &value_new,
                                std::vector<eigen_triplet_t> &basisvectors_triplets,
                                std::vector<double> &sqnorm_list) {
    auto state_iter = states.get<1>().find(state);

    size_t row_new;
    if (state_iter != states.get<1>().end()) {
        row_new = state_iter->idx;
    } else {
        row_new = states.size();
        states.emplace_back(row_new, state);
    }

    basisvectors_triplets.emplace_back(row_new, col_new, value_new);
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
