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

#include "SystemOne.hpp"
#include "Constants.hpp"
#include "MatrixElementCache.hpp"
#include "Symmetry.hpp"

#include <cctype>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

template <typename Scalar>
SystemOne<Scalar>::SystemOne(std::string species, MatrixElementCache &cache)
    : SystemBase<Scalar, StateOne>(cache), efield({{0, 0, 0}}), bfield({{0, 0, 0}}),
      diamagnetism(true), charge(0), ordermax(0), distance(std::numeric_limits<double>::max()),
      species(std::move(species)), sym_reflection(NA), sym_rotation({static_cast<float>(ARB)}) {}

template <typename Scalar>
SystemOne<Scalar>::SystemOne(std::string species, MatrixElementCache &cache, bool memory_saving)
    : SystemBase<Scalar, StateOne>(cache, memory_saving), efield({{0, 0, 0}}), bfield({{0, 0, 0}}),
      diamagnetism(true), charge(0), ordermax(0), distance(std::numeric_limits<double>::max()),
      species(std::move(species)), sym_reflection(NA), sym_rotation({static_cast<float>(ARB)}) {}

template <typename Scalar>
const std::string &SystemOne<Scalar>::getSpecies() const {
    return species;
}

template <typename Scalar>
void SystemOne<Scalar>::setEfield(std::array<double, 3> field) {
    this->onParameterChange();
    efield = field;

    // Transform the electric field into spherical coordinates
    this->changeToSphericalbasis(efield, efield_spherical);
}

template <typename Scalar>
void SystemOne<Scalar>::setBfield(std::array<double, 3> field) {
    this->onParameterChange();
    bfield = field;

    // Transform the magnetic field into spherical coordinates
    this->changeToSphericalbasis(bfield, bfield_spherical);

    diamagnetism_terms[{{0, +0}}] = bfield_spherical[+0] * bfield_spherical[+0] -
        bfield_spherical[+1] * bfield_spherical[-1] * 2.;
    diamagnetism_terms[{{2, +0}}] =
        bfield_spherical[+0] * bfield_spherical[+0] + bfield_spherical[+1] * bfield_spherical[-1];
    diamagnetism_terms[{{2, +1}}] = bfield_spherical[+0] * bfield_spherical[-1];
    diamagnetism_terms[{{2, -1}}] = bfield_spherical[+0] * bfield_spherical[+1];
    diamagnetism_terms[{{2, +2}}] = bfield_spherical[-1] * bfield_spherical[-1];
    diamagnetism_terms[{{2, -2}}] = bfield_spherical[+1] * bfield_spherical[+1];
}

template <typename Scalar>
void SystemOne<Scalar>::setEfield(std::array<double, 3> field, std::array<double, 3> to_z_axis,
                                  std::array<double, 3> to_y_axis) {
    this->rotateVector(field, to_z_axis, to_y_axis);
    this->setEfield(field);
}

template <typename Scalar>
void SystemOne<Scalar>::setBfield(std::array<double, 3> field, std::array<double, 3> to_z_axis,
                                  std::array<double, 3> to_y_axis) {
    this->rotateVector(field, to_z_axis, to_y_axis);
    this->setBfield(field);
}

template <typename Scalar>
void SystemOne<Scalar>::setEfield(std::array<double, 3> field, double alpha, double beta,
                                  double gamma) {
    this->rotateVector(field, alpha, beta, gamma);
    this->setEfield(field);
}

template <typename Scalar>
void SystemOne<Scalar>::setBfield(std::array<double, 3> field, double alpha, double beta,
                                  double gamma) {
    this->rotateVector(field, alpha, beta, gamma);
    this->setBfield(field);
}

template <typename Scalar>
void SystemOne<Scalar>::enableDiamagnetism(bool enable) {
    this->onParameterChange();
    diamagnetism = enable;
}
template <typename Scalar>
void SystemOne<Scalar>::setIonCharge(int c) {
    this->onParameterChange();
    charge = c;
}

template <typename Scalar>
void SystemOne<Scalar>::setRydIonOrder(unsigned int o) {
    this->onParameterChange();
    ordermax = o;
}

template <typename Scalar>
void SystemOne<Scalar>::setRydIonDistance(double d) {
    this->onParameterChange();
    distance = d;
}

template <typename Scalar>
void SystemOne<Scalar>::setConservedParityUnderReflection(parity_t parity) {
    this->onSymmetryChange();
    sym_reflection = parity;
    if (!this->isRefelectionAndRotationCompatible()) {
        throw std::runtime_error("The conserved parity under reflection is not compatible to the "
                                 "previously specified conserved momenta.");
    }
}

template <typename Scalar>
void SystemOne<Scalar>::setConservedMomentaUnderRotation(const std::set<float> &momenta) {
    if (momenta.count(static_cast<float>(ARB)) != 0 && momenta.size() > 1) {
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

template <typename Scalar>
void SystemOne<Scalar>::initializeBasis() {
    // If the basis is infinite, throw an error
    if (this->range_n.empty() &&
        (this->energy_min == std::numeric_limits<double>::lowest() ||
         this->energy_max == std::numeric_limits<double>::max())) {
        throw std::runtime_error(
            "The number of basis elements is infinite. The basis has to be restricted.");
    }

    ////////////////////////////////////////////////////////////////////
    /// Build one atom states //////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // TODO check whether symmetries are applicable
    // TODO check whether range_j, range_m is half-integer or integer valued

    float s = 0.5;
    if (std::isdigit(species.back()) != 0) {
        s = ((species.back() - '0') - 1) / 2.; // TODO think of a better solution
    }

    size_t idx = 0;
    std::vector<Eigen::Triplet<Scalar>>
        basisvectors_triplets; // TODO reserve states, basisvectors_triplets,
                               // hamiltonian_triplets

    std::vector<Eigen::Triplet<Scalar>> hamiltonian_triplets;

    /// Loop over specified quantum numbers ////////////////////////////

    std::set<int> range_adapted_n, range_adapted_l;
    std::set<float> range_adapted_j, range_adapted_m;

    if (this->range_n.empty()) {
        throw std::runtime_error(
            "The calculation of range_n via energy restrictions is not yet implemented."); // TODO
    }
    range_adapted_n = this->range_n;

    for (auto n : range_adapted_n) {

        if (this->range_l.empty()) {
            this->range(range_adapted_l, 0, n - 1);
        } else {
            range_adapted_l = this->range_l;
        }
        for (auto l : range_adapted_l) {
            if (l > n - 1 || l < 0) {
                continue;
            }

            if (this->range_j.empty()) {
                this->range(range_adapted_j, std::fabs(l - s), l + s);
            } else {
                range_adapted_j = this->range_j;
            }
            for (auto j : range_adapted_j) {
                if (std::fabs(j - l) > s || j < 0) {
                    continue;
                }

                double energy = StateOne(species, n, l, j, s).getEnergy(this->cache);
                if (!this->checkIsEnergyValid(energy)) {
                    continue;
                }

                if (this->range_m.empty()) {
                    this->range(range_adapted_m, -j, j);
                } else {
                    range_adapted_m = this->range_m;
                }

                // Consider rotation symmetry
                std::set<float> range_allowed_m;
                if (sym_rotation.count(static_cast<float>(ARB)) == 0) {
                    std::set_intersection(sym_rotation.begin(), sym_rotation.end(),
                                          range_adapted_m.begin(), range_adapted_m.end(),
                                          std::inserter(range_allowed_m, range_allowed_m.begin()));
                } else {
                    range_allowed_m = range_adapted_m;
                }

                for (auto m : range_allowed_m) {
                    if (std::fabs(m) > j) {
                        continue;
                    }

                    // Create state
                    StateOne state(species, n, l, j, m);

                    // Check whether reflection symmetry can be realized with the states available
                    if (sym_reflection != NA && state.getM() != 0 &&
                        range_allowed_m.count(-state.getM()) == 0) {
                        throw std::runtime_error("The momentum " + std::to_string(-state.getM()) +
                                                 " required by symmetries cannot be found.");
                    }

                    // Add symmetrized basis vectors
                    this->addSymmetrizedBasisvectors(state, idx, energy, basisvectors_triplets,
                                                     hamiltonian_triplets, sym_reflection);
                }
            }
        }
    }

    /// Loop over user-defined states //////////////////////////////////

    // Check that the user-defined states are not already contained in the list of states
    for (const auto &state : this->states_to_add) {
        if (this->states.template get<1>().find(state) != this->states.template get<1>().end()) {
            throw std::runtime_error("The state " + state.str() +
                                     " is already contained in the list of states.");
        }
        if (!state.isArtificial() && state.getSpecies() != species) {
            throw std::runtime_error("The state " + state.str() + " is of the wrong species.");
        }
    }

    // Add user-defined states
    for (const auto &state : this->states_to_add) {
        // Get energy of the state
        double energy = state.isArtificial() ? 0 : state.getEnergy(this->cache);

        // In case of artificial states, symmetries won't work
        auto sym_reflection_local = sym_reflection;
        auto sym_rotation_local = sym_rotation;
        if (state.isArtificial()) {
            if (sym_reflection_local != NA ||
                sym_rotation_local.count(static_cast<float>(ARB)) == 0) {
                std::cerr
                    << "WARNING: Only permutation symmetry can be applied to artificial states."
                    << std::endl;
            }
            sym_reflection_local = NA;
            sym_rotation_local = std::set<float>({static_cast<float>(ARB)});
        }

        // Consider rotation symmetry
        if (sym_rotation_local.count(static_cast<float>(ARB)) == 0 &&
            sym_rotation_local.count(state.getM()) == 0) {
            continue;
        }

        // Check whether reflection symmetry can be realized with the states available
        if (sym_reflection_local != NA && state.getM() != 0) {
            auto state_reflected = state.getReflected();
            if (this->states_to_add.find(state_reflected) == this->states_to_add.end()) {
                throw std::runtime_error("The state " + state_reflected.str() +
                                         " required by symmetries cannot be found.");
            }
        }

        // Add symmetrized basis vectors
        this->addSymmetrizedBasisvectors(state, idx, energy, basisvectors_triplets,
                                         hamiltonian_triplets, sym_reflection_local);
    }

    /// Build data /////////////////////////////////////////////////////

    this->basisvectors.resize(this->states.size(), idx);
    this->basisvectors.setFromTriplets(basisvectors_triplets.begin(), basisvectors_triplets.end());
    basisvectors_triplets.clear();

    this->hamiltonian.resize(idx, idx);
    this->hamiltonian.setFromTriplets(hamiltonian_triplets.begin(), hamiltonian_triplets.end());
    hamiltonian_triplets.clear();
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to calculate the interaction /////
////////////////////////////////////////////////////////////////////

template <typename Scalar>
void SystemOne<Scalar>::initializeInteraction() {
    ////////////////////////////////////////////////////////////////////
    /// Prepare the calculation of the interaction /////////////////////
    ////////////////////////////////////////////////////////////////////

    // Check if something to do
    double tolerance = 1e-24;

    std::vector<int> erange, brange;
    std::vector<std::array<int, 2>> drange;
    std::vector<int> orange;
    for (const auto &entry : efield_spherical) {
        if (entry.first < 0) {
            continue;
        }
        if (std::abs(entry.second) > tolerance &&
            interaction_efield.find(-entry.first) == interaction_efield.end()) {
            erange.push_back(entry.first);
        }
    }
    for (const auto &entry : bfield_spherical) {
        if (entry.first < 0) {
            continue;
        }
        if (std::abs(entry.second) > tolerance &&
            interaction_bfield.find(-entry.first) == interaction_bfield.end()) {
            brange.push_back(entry.first);
        }
    }
    for (const auto &entry : diamagnetism_terms) {
        if (entry.first[1] < 0) {
            continue;
        }
        if (diamagnetism && std::abs(entry.second) > tolerance &&
            interaction_diamagnetism.find(entry.first) == interaction_diamagnetism.end()) {
            drange.push_back(entry.first);
        }
    }

    if (charge != 0) {
        for (unsigned int order = 1; order <= ordermax; ++order) {
            if (interaction_multipole.find(order) == interaction_multipole.end()) {
                orange.push_back(order);
            }
        }
    }
    // Return if there is nothing to do
    if (erange.empty() && brange.empty() && drange.empty() && orange.empty()) {
        return;
    }

    // Precalculate matrix elements
    auto states_converted = this->getStates();
    for (const auto &i : erange) {
        this->cache.precalculateElectricMomentum(states_converted, i);
        if (i != 0) {
            this->cache.precalculateElectricMomentum(states_converted, -i);
        }
    }
    for (const auto &i : brange) {
        this->cache.precalculateMagneticMomentum(states_converted, i);
        if (i != 0) {
            this->cache.precalculateMagneticMomentum(states_converted, -i);
        }
    }
    for (const auto &i : drange) {
        this->cache.precalculateDiamagnetism(states_converted, i[0], i[1]);
        if (i[1] != 0) {
            this->cache.precalculateDiamagnetism(states_converted, i[0], -i[1]);
        }
    }
    if (charge != 0) {
        for (unsigned int order = 1; order <= ordermax; ++order) {
            this->cache.precalculateMultipole(states_converted, order);
        }
    }

    ////////////////////////////////////////////////////////////////////
    /// Calculate the interaction in the canonical basis ///////////////
    ////////////////////////////////////////////////////////////////////

    std::unordered_map<int, std::vector<Eigen::Triplet<Scalar>>>
        interaction_efield_triplets; // TODO reserve
    std::unordered_map<int, std::vector<Eigen::Triplet<Scalar>>>
        interaction_bfield_triplets; // TODO reserve
    std::unordered_map<std::array<int, 2>, std::vector<Eigen::Triplet<Scalar>>,
                       utils::hash<std::array<int, 2>>>
        interaction_diamagnetism_triplets; // TODO reserve
    std::unordered_map<int, std::vector<Eigen::Triplet<Scalar>>>
        interaction_multipole_triplets; // TODO reserve
    // Loop over column entries
    for (const auto &c : this->states) { // TODO parallelization
        if (c.state.isArtificial()) {
            continue;
        }

        // Loop over row entries
        for (const auto &r : this->states) {
            if (r.state.isArtificial()) {
                continue;
            }

            // E-field interaction
            for (const auto &i : erange) {
                if (i == 0 && r.idx < c.idx) {
                    continue;
                }

                if (selectionRulesMultipoleNew(r.state, c.state, 1, i)) {
                    Scalar value = this->cache.getElectricDipole(r.state, c.state);
                    this->addTriplet(interaction_efield_triplets[i], r.idx, c.idx, value);
                    break; // because for the other operators, the selection rule for the magnetic
                           // quantum numbers will not be fulfilled
                }
            }

            // B-field interaction
            for (const auto &i : brange) {
                if (i == 0 && r.idx < c.idx) {
                    continue;
                }

                if (selectionRulesMomentumNew(r.state, c.state, i)) {
                    Scalar value = this->cache.getMagneticDipole(r.state, c.state);
                    this->addTriplet(interaction_bfield_triplets[i], r.idx, c.idx, value);
                    break; // because for the other operators, the selection rule for the magnetic
                           // quantum numbers will not be fulfilled
                }
            }

            // Diamagnetic interaction
            for (const auto &i : drange) {
                if (i[1] == 0 && r.idx < c.idx) {
                    continue;
                }

                if (selectionRulesMultipoleNew(r.state, c.state, i[0], i[1])) {
                    Scalar value = 1. / (8 * electron_rest_mass) *
                        this->cache.getDiamagnetism(r.state, c.state, i[0]);
                    this->addTriplet(interaction_diamagnetism_triplets[i], r.idx, c.idx, value);
                }
            }

            // Multipole interaction with an ion
            if (charge != 0) {
                int q = r.state.getM() - c.state.getM();
                if (q == 0) { // total momentum consreved
                    for (const auto &order : orange) {
                        if (selectionRulesMultipoleNew(r.state, c.state, order)) {
                            double val = -coulombs_constant * elementary_charge *
                                this->cache.getElectricMultipole(r.state, c.state, order);
                            this->addTriplet(interaction_multipole_triplets[order], r.idx, c.idx,
                                             val);
                        }
                    }
                }
            }
        }
    }
    ////////////////////////////////////////////////////////////////////
    /// Build and transform the interaction to the used basis //////////
    ////////////////////////////////////////////////////////////////////

    for (const auto &i : erange) {
        interaction_efield[i].resize(this->states.size(), this->states.size());
        interaction_efield[i].setFromTriplets(interaction_efield_triplets[i].begin(),
                                              interaction_efield_triplets[i].end());
        interaction_efield_triplets[i].clear();

        if (i == 0) {
            interaction_efield[i] = this->basisvectors.adjoint() *
                interaction_efield[i].template selfadjointView<Eigen::Lower>() * this->basisvectors;
        } else {
            interaction_efield[i] =
                this->basisvectors.adjoint() * interaction_efield[i] * this->basisvectors;
            interaction_efield[-i] = std::pow(-1, i) * interaction_efield[i].adjoint();
        }
    }

    for (const auto &i : brange) {
        interaction_bfield[i].resize(this->states.size(), this->states.size());
        interaction_bfield[i].setFromTriplets(interaction_bfield_triplets[i].begin(),
                                              interaction_bfield_triplets[i].end());
        interaction_bfield_triplets[i].clear();

        if (i == 0) {
            interaction_bfield[i] = this->basisvectors.adjoint() *
                interaction_bfield[i].template selfadjointView<Eigen::Lower>() * this->basisvectors;
        } else {
            interaction_bfield[i] =
                this->basisvectors.adjoint() * interaction_bfield[i] * this->basisvectors;
            interaction_bfield[-i] = std::pow(-1, i) * interaction_bfield[i].adjoint();
        }
    }

    for (const auto &i : drange) {
        interaction_diamagnetism[i].resize(this->states.size(), this->states.size());
        interaction_diamagnetism[i].setFromTriplets(interaction_diamagnetism_triplets[i].begin(),
                                                    interaction_diamagnetism_triplets[i].end());
        interaction_diamagnetism_triplets[i].clear();

        if (i[1] == 0) {
            interaction_diamagnetism[i] = this->basisvectors.adjoint() *
                interaction_diamagnetism[i].template selfadjointView<Eigen::Lower>() *
                this->basisvectors;
        } else {
            interaction_diamagnetism[i] =
                this->basisvectors.adjoint() * interaction_diamagnetism[i] * this->basisvectors;
            interaction_diamagnetism[{{i[0], -i[1]}}] =
                std::pow(-1, i[1]) * interaction_diamagnetism[i].adjoint();
        }
    }
    if (charge != 0) {
        for (const auto &i : orange) {
            interaction_multipole[i].resize(this->states.size(), this->states.size());
            interaction_multipole[i].setFromTriplets(interaction_multipole_triplets[i].begin(),
                                                     interaction_multipole_triplets[i].end());
            interaction_multipole_triplets[i].clear();
            if (i == 0) {
                interaction_multipole[i] = this->basisvectors.adjoint() *
                    interaction_multipole[i].template selfadjointView<Eigen::Lower>() *
                    this->basisvectors;
            } else {
                interaction_multipole[i] =
                    this->basisvectors.adjoint() * interaction_multipole[i] * this->basisvectors;
                interaction_multipole[-i] = std::pow(-1, i) * interaction_multipole[i].adjoint();
            }
        }
    }
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to construct Hamiltonian /////////
////////////////////////////////////////////////////////////////////

template <typename Scalar>
void SystemOne<Scalar>::addInteraction() {
    // Build the total Hamiltonian
    double tolerance = 1e-24;

    if (std::abs(efield_spherical[+0]) > tolerance) {
        this->hamiltonian -= interaction_efield[+0] * efield_spherical[+0];
    }
    if (std::abs(efield_spherical[-1]) > tolerance) {
        this->hamiltonian += interaction_efield[+1] * efield_spherical[-1];
    }
    if (std::abs(efield_spherical[+1]) > tolerance) {
        this->hamiltonian += interaction_efield[-1] * efield_spherical[+1];
    }
    if (std::abs(bfield_spherical[+0]) > tolerance) {
        this->hamiltonian -= interaction_bfield[+0] * bfield_spherical[+0];
    }
    if (std::abs(bfield_spherical[-1]) > tolerance) {
        this->hamiltonian += interaction_bfield[+1] * bfield_spherical[-1];
    }
    if (std::abs(bfield_spherical[+1]) > tolerance) {
        this->hamiltonian += interaction_bfield[-1] * bfield_spherical[+1];
    }

    if (diamagnetism && std::abs(diamagnetism_terms[{{0, +0}}]) > tolerance) {
        this->hamiltonian += interaction_diamagnetism[{{0, +0}}] * diamagnetism_terms[{{0, +0}}];
    }
    if (diamagnetism && std::abs(diamagnetism_terms[{{2, +0}}]) > tolerance) {
        this->hamiltonian -= interaction_diamagnetism[{{2, +0}}] * diamagnetism_terms[{{2, +0}}];
    }
    if (diamagnetism && std::abs(diamagnetism_terms[{{2, +1}}]) > tolerance) {
        this->hamiltonian +=
            interaction_diamagnetism[{{2, +1}}] * diamagnetism_terms[{{2, +1}}] * std::sqrt(3);
    }
    if (diamagnetism && std::abs(diamagnetism_terms[{{2, -1}}]) > tolerance) {
        this->hamiltonian +=
            interaction_diamagnetism[{{2, -1}}] * diamagnetism_terms[{{2, -1}}] * std::sqrt(3);
    }
    if (diamagnetism && std::abs(diamagnetism_terms[{{2, +2}}]) > tolerance) {
        this->hamiltonian -=
            interaction_diamagnetism[{{2, +2}}] * diamagnetism_terms[{{2, +2}}] * std::sqrt(1.5);
    }
    if (diamagnetism && std::abs(diamagnetism_terms[{{2, -2}}]) > tolerance) {
        this->hamiltonian -=
            interaction_diamagnetism[{{2, -2}}] * diamagnetism_terms[{{2, -2}}] * std::sqrt(1.5);
    }
    if (charge != 0 && distance != std::numeric_limits<double>::max()) {
        for (unsigned int order = 1; order <= ordermax; ++order) {
            double powerlaw = 1. / std::pow(distance, order + 1);
            this->hamiltonian += interaction_multipole[order] * charge * powerlaw;
        }
    }
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to transform the interaction /////
////////////////////////////////////////////////////////////////////

template <typename Scalar>
void SystemOne<Scalar>::transformInteraction(const Eigen::SparseMatrix<Scalar> &transformator) {
    for (auto &entry : interaction_efield) {
        entry.second = transformator.adjoint() * entry.second * transformator;
    }
    for (auto &entry : interaction_bfield) {
        entry.second = transformator.adjoint() * entry.second * transformator;
    }
    for (auto &entry : interaction_diamagnetism) {
        entry.second = transformator.adjoint() * entry.second * transformator; // NOLINT
    }
    for (auto &entry : interaction_multipole) {
        entry.second = transformator.adjoint() * entry.second * transformator; // NOLINT
    }
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to delete the interaction ////////
////////////////////////////////////////////////////////////////////

template <typename Scalar>
void SystemOne<Scalar>::deleteInteraction() {
    interaction_efield.clear();
    interaction_bfield.clear();
    interaction_diamagnetism.clear();
    interaction_multipole.clear();
}

////////////////////////////////////////////////////////////////////
/// Methods that allows base class to rotate states ////////////////
////////////////////////////////////////////////////////////////////

template <typename Scalar>
Eigen::SparseMatrix<Scalar>
SystemOne<Scalar>::rotateStates(const std::vector<size_t> &states_indices, double alpha,
                                double beta, double gamma) {
    // Initialize Wigner D matrix
    WignerD wigner;

    // Rotate state
    std::vector<Eigen::Triplet<Scalar>> states_rotated_triplets;
    states_rotated_triplets.reserve(
        std::min(static_cast<size_t>(10), this->states.size()) *
        states_indices.size()); // TODO std::min( 2*jmax+1, states.size() ) * states_indices.size()

    size_t current = 0;
    for (auto const &idx : states_indices) {
        this->addRotated(this->states[idx].state, current++, states_rotated_triplets, wigner, alpha,
                         beta, gamma);
    }

    Eigen::SparseMatrix<Scalar> states_rotated(this->states.size(), states_indices.size());
    states_rotated.setFromTriplets(states_rotated_triplets.begin(), states_rotated_triplets.end());
    states_rotated_triplets.clear();

    return states_rotated;
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar> SystemOne<Scalar>::buildStaterotator(double alpha, double beta,
                                                                 double gamma) {
    // Initialize Wigner D matrix
    WignerD wigner;

    // Build rotator
    std::vector<Eigen::Triplet<Scalar>> rotator_triplets;
    rotator_triplets.reserve(
        std::min(static_cast<size_t>(10), this->states.size()) *
        this->states.size()); // TODO std::min( 2*jmax+1, states.size() ) * states.size()

    for (auto const &entry : this->states) {
        this->addRotated(entry.state, entry.idx, rotator_triplets, wigner, alpha, beta, gamma);
    }

    Eigen::SparseMatrix<Scalar> rotator(this->states.size(), this->states.size());
    rotator.setFromTriplets(rotator_triplets.begin(), rotator_triplets.end()); // NOLINT
    rotator_triplets.clear();

    return rotator;
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to combine systems ///////////////
////////////////////////////////////////////////////////////////////

template <typename Scalar>
void SystemOne<Scalar>::incorporate(SystemBase<Scalar, StateOne> &system) {
    // Combine parameters
    if (species != dynamic_cast<SystemOne<Scalar> &>(system).species) {
        throw std::runtime_error(
            "The value of the variable 'element' must be the same for both systems.");
    }
    if (efield != dynamic_cast<SystemOne<Scalar> &>(system).efield) {
        throw std::runtime_error("The value of the variable 'distance' must be the same for both "
                                 "systems."); // implies that efield_spherical is the same, too
    }
    if (bfield != dynamic_cast<SystemOne<Scalar> &>(system).bfield) {
        throw std::runtime_error("The value of the variable 'angle' must be the same for both "
                                 "systems."); // implies that
                                              // bfield_spherical
                                              // is the same,
                                              // too
    }
    if (diamagnetism != dynamic_cast<SystemOne<Scalar> &>(system).diamagnetism) {
        throw std::runtime_error(
            "The value of the variable 'ordermax' must be the same for both systems.");
    }

    // Combine symmetries
    unsigned int num_different_symmetries = 0;
    if (sym_reflection != dynamic_cast<SystemOne<Scalar> &>(system).sym_reflection) {
        sym_reflection = NA;
        ++num_different_symmetries;
    }
    if (!(sym_rotation.size() == dynamic_cast<SystemOne<Scalar> &>(system).sym_rotation.size() &&
          std::equal(sym_rotation.begin(), sym_rotation.end(),
                     dynamic_cast<SystemOne<Scalar> &>(system).sym_rotation.begin()))) {
        if (sym_rotation.count(static_cast<float>(ARB)) != 0 ||
            dynamic_cast<SystemOne<Scalar> &>(system).sym_rotation.count(static_cast<float>(ARB)) !=
                0) {
            sym_rotation = {static_cast<float>(ARB)};
        } else {
            sym_rotation.insert(dynamic_cast<SystemOne<Scalar> &>(system).sym_rotation.begin(),
                                dynamic_cast<SystemOne<Scalar> &>(system).sym_rotation.end());
        }
        ++num_different_symmetries;
    }
    if (num_different_symmetries > 1) {
        std::cerr << "Warning: The systems differ in more than one symmetry. For the combined "
                     "system, the notion of symmetries might be meaningless."
                  << std::endl;
    }

    // Clear cached interaction
    this->deleteInteraction();
}

////////////////////////////////////////////////////////////////////
/// Utility methods ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

template <typename Scalar>
void SystemOne<Scalar>::addSymmetrizedBasisvectors(
    const StateOne &state, size_t &idx, const double &energy,
    std::vector<Eigen::Triplet<Scalar>> &basisvectors_triplets,
    std::vector<Eigen::Triplet<Scalar>> &hamiltonian_triplets, parity_t &sym_reflection_local) {
    // In case of reflection symmetry, skip half of the basis vectors
    if (sym_reflection_local != NA && state.getM() != 0) {
        if (state.getM() < 0) {
            return;
        }
    }

    // Store the energy of the unperturbed one atom state
    hamiltonian_triplets.emplace_back(idx, idx, energy);

    // Adapt the normalization if required by symmetries
    Scalar value = 1;
    if (sym_reflection_local != NA && state.getM() != 0) {
        value /= std::sqrt(2);
    }

    // Add an entry to the current basis vector
    this->addBasisvectors(state, idx, value, basisvectors_triplets);

    // Add further entries to the current basis vector if required by symmetries
    if (sym_reflection_local != NA && state.getM() != 0) {
        value *= std::pow(-1, state.getL() + state.getM() - state.getJ()) *
            utils::imaginary_unit<Scalar>();
        value *= (sym_reflection_local == EVEN) ? 1 : -1;
        // S_y is invariant under reflection through xz-plane
        // TODO is the s quantum number of importance here?
        this->addBasisvectors(state.getReflected(), idx, value, basisvectors_triplets);
    }

    ++idx;
}

template <typename Scalar>
void SystemOne<Scalar>::addBasisvectors(
    const StateOne &state, const size_t &idx, const Scalar &value,
    std::vector<Eigen::Triplet<Scalar>> &basisvectors_triplets) {
    auto state_iter = this->states.template get<1>().find(state);

    size_t row;
    if (state_iter != this->states.template get<1>().end()) {
        row = state_iter->idx;
    } else {
        row = this->states.size();
        this->states.push_back(enumerated_state<StateOne>(row, state));
    }

    basisvectors_triplets.emplace_back(row, idx, value);
}

template <typename Scalar>
void SystemOne<Scalar>::changeToSphericalbasis(std::array<double, 3> field,
                                               std::unordered_map<int, double> &field_spherical) {
    if (field[1] != 0) {
        throw std::runtime_error(
            "For fields with non-zero y-coordinates, a complex data type is needed.");
    }
    field_spherical[1] = -field[0] / std::sqrt(2);
    field_spherical[-1] = field[0] / std::sqrt(2);
    field_spherical[0] = field[2];
}

template <typename Scalar>
void SystemOne<Scalar>::changeToSphericalbasis(
    std::array<double, 3> field, std::unordered_map<int, std::complex<double>> &field_spherical) {
    field_spherical[1] = std::complex<double>(-field[0] / std::sqrt(2), -field[1] / std::sqrt(2));
    field_spherical[-1] = std::complex<double>(field[0] / std::sqrt(2), -field[1] / std::sqrt(2));
    field_spherical[0] = std::complex<double>(field[2], 0);
}

template <typename Scalar>
void SystemOne<Scalar>::addTriplet(std::vector<Eigen::Triplet<Scalar>> &triplets,
                                   const size_t r_idx, const size_t c_idx, const Scalar val) {
    triplets.emplace_back(r_idx, c_idx, val);
}

template <typename Scalar>
void SystemOne<Scalar>::rotateVector(std::array<double, 3> &field, std::array<double, 3> &to_z_axis,
                                     std::array<double, 3> &to_y_axis) {
    auto field_mapped = Eigen::Map<Eigen::Matrix<double, 3, 1>>(&field[0]);

    if (field_mapped.norm() != 0) {
        Eigen::Matrix<double, 3, 3> rotator = this->buildRotator(to_z_axis, to_y_axis);
        field_mapped = rotator.transpose() * field_mapped;
    }
}

template <typename Scalar>
void SystemOne<Scalar>::rotateVector(std::array<double, 3> &field, double alpha, double beta,
                                     double gamma) {
    auto field_mapped = Eigen::Map<Eigen::Matrix<double, 3, 1>>(&field[0]);

    if (field_mapped.norm() != 0) {
        Eigen::Matrix<double, 3, 3> rotator = this->buildRotator(alpha, beta, gamma);
        field_mapped = rotator.transpose() * field_mapped;
    }
}

template <typename Scalar>
bool SystemOne<Scalar>::isRefelectionAndRotationCompatible() {
    if (sym_rotation.count(static_cast<float>(ARB)) != 0 || sym_reflection == NA) {
        return true;
    }

    for (const auto &s : sym_rotation) {
        if (sym_rotation.count(-s) == 0) {
            return false;
        }
    }

    return true;
}

template class SystemOne<std::complex<double>>;
template class SystemOne<double>;
