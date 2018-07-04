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

#include <cctype>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include <unordered_set>
#include <type_traits>

SystemOne::SystemOne(std::string  species, MatrixElementCache &cache)
    : SystemBase(cache), efield({{0,0,0}}), bfield({{0,0,0}}), diamagnetism(false), species(std::move(species)), sym_reflection(NA), sym_rotation({static_cast<float>(ARB)}) {
}

SystemOne::SystemOne(std::string  species, MatrixElementCache &cache, bool memory_saving)
    : SystemBase(cache, memory_saving), efield({{0,0,0}}), bfield({{0,0,0}}), diamagnetism(false), species(std::move(species)), sym_reflection(NA), sym_rotation({static_cast<float>(ARB)}) {
}

const std::string& SystemOne::getElement() const {
    return species;
}

void SystemOne::setEfield(std::array<double, 3> field) {
    this->onParameterChange();
    efield = field;

    // Transform the electric field into spherical coordinates
    this->changeToSphericalbasis(efield, efield_spherical);
}

void SystemOne::setBfield(std::array<double, 3> field) {
    this->onParameterChange();
    bfield = field;

    // Transform the magnetic field into spherical coordinates
    this->changeToSphericalbasis(bfield, bfield_spherical);

    diamagnetism_terms[{{0,+0}}] = bfield_spherical[+0]*bfield_spherical[+0]-bfield_spherical[+1]*bfield_spherical[-1]*2.;
    diamagnetism_terms[{{2,+0}}] = bfield_spherical[+0]*bfield_spherical[+0]+bfield_spherical[+1]*bfield_spherical[-1];
    diamagnetism_terms[{{2,+1}}] = bfield_spherical[+0]*bfield_spherical[-1];
    diamagnetism_terms[{{2,-1}}] = bfield_spherical[+0]*bfield_spherical[+1];
    diamagnetism_terms[{{2,+2}}] = bfield_spherical[-1]*bfield_spherical[-1];
    diamagnetism_terms[{{2,-2}}] = bfield_spherical[+1]*bfield_spherical[+1];
}

void SystemOne::setEfield(std::array<double, 3> field, std::array<double, 3> to_z_axis, std::array<double, 3> to_y_axis) {
    this->rotateVector(field, to_z_axis, to_y_axis);
    this->setEfield(field);
}

void SystemOne::setBfield(std::array<double, 3> field, std::array<double, 3> to_z_axis, std::array<double, 3> to_y_axis) {
    this->rotateVector(field, to_z_axis, to_y_axis);
    this->setBfield(field);
}


void SystemOne::setEfield(std::array<double, 3> field, double alpha, double beta, double gamma) {
    this->rotateVector(field, alpha, beta, gamma);
    this->setEfield(field);
}

void SystemOne::setBfield(std::array<double, 3> field, double alpha, double beta, double gamma) {
    this->rotateVector(field, alpha, beta, gamma);
    this->setBfield(field);
}

void SystemOne::setDiamagnetism(bool enable) {
    this->onParameterChange();
    diamagnetism = enable;
}

void SystemOne::setConservedParityUnderReflection(parity_t parity) {
    this->onSymmetryChange();
    sym_reflection = parity;
    if (!this->isRefelectionAndRotationCompatible()) { throw std::runtime_error("The conserved parity under reflection is not compatible to the previously specified conserved momenta.");
}
}

void SystemOne::setConservedMomentaUnderRotation(const std::set<float>& momenta) {
    if (momenta.count(static_cast<float>(ARB))!=0 && momenta.size() > 1) { throw std::runtime_error("If ARB (=arbitrary momentum) is specified, momenta must not be passed explicitely.");
}
    this->onSymmetryChange();
    sym_rotation = momenta;
    if (!this->isRefelectionAndRotationCompatible()) { throw std::runtime_error("The conserved momenta are not compatible to the previously specified conserved parity under reflection.");
}
}

////////////////////////////////////////////////////////////////////
/// Explicit template specializations (has to be at the beginning) /
////////////////////////////////////////////////////////////////////

template<>
double SystemOne::imaginaryUnit() {
    throw std::runtime_error( "For operations that invoke the imaginary number, a complex data type is needed." );
}

template<>
std::complex<double> SystemOne::imaginaryUnit() {
    return {0,1};
}

template<>
double SystemOne::convert(const std::complex<double> &val) {
    return val.real();
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to initialize Basis //////////////
////////////////////////////////////////////////////////////////////

void SystemOne::initializeBasis()
{
    // If the basis is infinite, throw an error
    if (range_n.empty() && (energy_min == std::numeric_limits<double>::lowest() || energy_max == std::numeric_limits<double>::max())) {
        throw std::runtime_error( "The number of basis elements is infinite. The basis has to be restricted." );
    }

    ////////////////////////////////////////////////////////////////////
    /// Build one atom states //////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // TODO check whether symmetries are applicable
    // TODO check whether range_j, range_m is half-integer or integer valued

    float s = 0.5;
    if (std::isdigit(species.back()) != 0) { s = (std::atoi(&species.back())-1)/2.; // TODO think of a better solution
}

    size_t idx = 0;
    std::vector<eigen_triplet_t> coefficients_triplets; // TODO reserve states, coefficients_triplets, hamiltonianmatrix_triplets
    std::vector<eigen_triplet_t> hamiltonianmatrix_triplets;
    std::set<int> range_adapted_n, range_adapted_l;
    std::set<float> range_adapted_j, range_adapted_m;

    if (range_n.empty()) {
        throw std::runtime_error( "The calculation of range_n via energy restrictions is not yet implemented." ); // TODO
    } 
        range_adapted_n = range_n;
    
    for (auto n : range_adapted_n) {

        if (range_l.empty()) {
            this->range(range_adapted_l, 0, n-1);
        } else {
            range_adapted_l = range_l;
        }
        for (auto l : range_adapted_l) {
            if (l > n-1 || l < 0) { continue;
}

            if (range_j.empty()) {
                this->range(range_adapted_j, std::fabs(l-s), l+s);
            } else {
                range_adapted_j = range_j;
            }
            for (auto j : range_adapted_j) {
                if (std::fabs(j-l) > s || j < 0) { continue;
}

                double energy = StateOne(species,n,l,j,s).getEnergy();
                if (!checkIsEnergyValid(energy)) { continue;
}

                if (range_m.empty()) {
                    this->range(range_adapted_m, -j, j);
                } else {
                    range_adapted_m = range_m;
                }

                // Consider rotation symmetry
                std::set<float> range_allowed_m;
                if (sym_rotation.count(static_cast<float>(ARB)) == 0) {
                    std::set_intersection(sym_rotation.begin(), sym_rotation.end(), range_adapted_m.begin(), range_adapted_m.end(), std::inserter(range_allowed_m, range_allowed_m.begin()));
                } else {
                    range_allowed_m = range_adapted_m;
                }

                for (auto m : range_allowed_m) {
                    if (std::fabs(m) > j) { continue;
}

                    // In case of reflection symmetry, skip half of the basis vectors
                    if (sym_reflection != NA && m < 0) { continue;
}

                    // Store the energy of the unperturbed one atom state
                    hamiltonianmatrix_triplets.emplace_back(idx,idx,energy);

                    // Adapt the normalization if required by symmetries
                    scalar_t value = 1;
                    if (sym_reflection != NA) {
                        value /= std::sqrt(2);
                    }

                    // Add an entry to the current basis vector
                    this->addCoefficient(StateOne(species,n,l,j,m), idx, value, coefficients_triplets);

                    // Add further entries to the current basis vector if required by symmetries
                    if (sym_reflection != NA) {
                        value *= (sym_reflection == EVEN) ? std::pow(-1,l+m-j)*this->imaginaryUnit<scalar_t>() : -std::pow(-1,l+m-j)*this->imaginaryUnit<scalar_t>(); // S_y is invariant under reflection through xz-plane
                        this->addCoefficient(StateOne(species,n,l,j,-m), idx, value, coefficients_triplets);
                    }

                    ++idx;
                }
            }
        }
    }

    // Build data
    coefficients.resize(states.size(),idx);
    coefficients.setFromTriplets(coefficients_triplets.begin(), coefficients_triplets.end());
    coefficients_triplets.clear();

    hamiltonianmatrix.resize(idx,idx);
    hamiltonianmatrix.setFromTriplets(hamiltonianmatrix_triplets.begin(), hamiltonianmatrix_triplets.end());
    hamiltonianmatrix_triplets.clear();
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to calculate the interaction /////
////////////////////////////////////////////////////////////////////

void SystemOne::initializeInteraction() {
    ////////////////////////////////////////////////////////////////////
    /// Prepare the calculation of the interaction /////////////////////
    ////////////////////////////////////////////////////////////////////

    // Check if something to do
    double tolerance = 1e-24;

    std::vector<int> erange, brange;
    std::vector<std::array<int, 2>> drange;

    for (const auto &entry : efield_spherical) {
        if (std::abs(entry.second) > tolerance && interaction_efield.find(-entry.first) == interaction_efield.end()) {
            erange.push_back(-entry.first); // TODO use entry.first without a minus
        }
    }
    for (const auto &entry : bfield_spherical) {
        if (std::abs(entry.second) > tolerance && interaction_bfield.find(-entry.first) == interaction_bfield.end()) {
            brange.push_back(-entry.first); // TODO use entry.first without a minus
        }
    }
    for (const auto &entry : diamagnetism_terms) {
        if (diamagnetism && std::abs(entry.second) > tolerance && interaction_diamagnetism.find(entry.first) == interaction_diamagnetism.end()) {
            drange.push_back(entry.first);
        }
    }

    if (erange.empty() && brange.empty() && drange.empty()) { return;
}

    // Precalculate matrix elements
    auto states_converted = this->getStates(); // TODO remove
    for (const auto &i : erange) { cache.precalculateElectricMomentum(states_converted, i);
}
    for (const auto &i : brange) { cache.precalculateMagneticMomentum(states_converted, i);
}
    for (const auto &i : drange) { cache.precalculateDiamagnetism(states_converted, i[0], i[1]);
}

    ////////////////////////////////////////////////////////////////////
    /// Calculate the interaction in the canonical basis ///////////////
    ////////////////////////////////////////////////////////////////////

    std::unordered_map<int, std::vector<eigen_triplet_t>> interaction_efield_triplets; // TODO reserve
    std::unordered_map<int, std::vector<eigen_triplet_t>> interaction_bfield_triplets;// TODO reserve
    std::unordered_map<std::array<int, 2>, std::vector<eigen_triplet_t>> interaction_diamagnetism_triplets; // TODO reserve

    // Loop over column entries
    for (const auto &c: states) { // TODO parallelization

        if (c.state.species.empty()) { continue; // TODO artifical states TODO [dummystates]
}

        // Loop over row entries
        for (const auto &r: states) {

            if (r.state.species.empty()) { continue; // TODO artifical states TODO [dummystates]
}
            //if (r.idx < c.idx) continue; // TODO modify addTriplet so that skipping half of the entries work

            // E-field interaction
            for (const auto &i : erange) {
                if (selectionRulesMultipoleNew(r.state, c.state, 1, i)) {
                    scalar_t value = cache.getElectricDipole(r.state, c.state);
                    this->addTriplet(interaction_efield_triplets[i], r.idx, c.idx, value);
                    break; // because for the other operators, the selection rule for the magnetic quantum numbers will not be fulfilled
                }
            }

            // B-field interaction
            for (const auto &i : brange) {
                if (selectionRulesMomentumNew(r.state, c.state, i)) {
                    scalar_t value = cache.getMagneticDipole(r.state, c.state);
                    this->addTriplet(interaction_bfield_triplets[i], r.idx, c.idx, value);
                    break; // because for the other operators, the selection rule for the magnetic quantum numbers will not be fulfilled
                }
            }

            // Diamagnetic interaction
            for (const auto &i : drange) {
                if (selectionRulesMultipoleNew(r.state, c.state, i[0], i[1])) {
                    scalar_t value = 1./(8*electron_rest_mass) * cache.getDiamagnetism(r.state, c.state, i[0]);
                    this->addTriplet(interaction_diamagnetism_triplets[i], r.idx, c.idx, value);
                }
            }
        }
    }

    ////////////////////////////////////////////////////////////////////
    /// Build and transform the interaction to the used basis //////////
    ////////////////////////////////////////////////////////////////////

    for (const auto &i : erange) {
        interaction_efield[i].resize(states.size(),states.size());
        interaction_efield[i].setFromTriplets(interaction_efield_triplets[i].begin(), interaction_efield_triplets[i].end());
        interaction_efield_triplets[i].clear();

        interaction_efield[i] = coefficients.adjoint()*interaction_efield[i]*coefficients;
    }

    for (const auto &i : brange) {
        interaction_bfield[i].resize(states.size(),states.size());
        interaction_bfield[i].setFromTriplets(interaction_bfield_triplets[i].begin(), interaction_bfield_triplets[i].end());
        interaction_bfield_triplets[i].clear();

        interaction_bfield[i] = coefficients.adjoint()*interaction_bfield[i]*coefficients;
    }

    for (const auto &i : drange) {
        interaction_diamagnetism[i].resize(states.size(),states.size());
        interaction_diamagnetism[i].setFromTriplets(interaction_diamagnetism_triplets[i].begin(), interaction_diamagnetism_triplets[i].end());
        interaction_diamagnetism_triplets[i].clear();

        interaction_diamagnetism[i] = coefficients.adjoint()*interaction_diamagnetism[i]*coefficients;
    }
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to construct Hamiltonian /////////
////////////////////////////////////////////////////////////////////

void SystemOne::addInteraction() {
    // Build the total Hamiltonian
    double tolerance = 1e-24;

    if (std::abs(efield_spherical[+0]) > tolerance) { hamiltonianmatrix -= interaction_efield[+0]*efield_spherical[+0];
}
    if (std::abs(efield_spherical[-1]) > tolerance) { hamiltonianmatrix += interaction_efield[+1]*efield_spherical[-1];
}
    if (std::abs(efield_spherical[+1]) > tolerance) { hamiltonianmatrix += interaction_efield[-1]*efield_spherical[+1];
}
    if (std::abs(bfield_spherical[+0]) > tolerance) { hamiltonianmatrix -= interaction_bfield[+0]*bfield_spherical[+0];
}
    if (std::abs(bfield_spherical[-1]) > tolerance) { hamiltonianmatrix += interaction_bfield[+1]*bfield_spherical[-1];
}
    if (std::abs(bfield_spherical[+1]) > tolerance) { hamiltonianmatrix += interaction_bfield[-1]*bfield_spherical[+1];
}

    if (diamagnetism && std::abs(diamagnetism_terms[{{0,+0}}]) > tolerance) { hamiltonianmatrix += interaction_diamagnetism[{{0,+0}}]*diamagnetism_terms[{{0,+0}}];
}
    if (diamagnetism && std::abs(diamagnetism_terms[{{2,+0}}]) > tolerance) { hamiltonianmatrix -= interaction_diamagnetism[{{2,+0}}]*diamagnetism_terms[{{2,+0}}];
}
    if (diamagnetism && std::abs(diamagnetism_terms[{{2,+1}}]) > tolerance) { hamiltonianmatrix += interaction_diamagnetism[{{2,+1}}]*diamagnetism_terms[{{2,+1}}]*std::sqrt(3);
}
    if (diamagnetism && std::abs(diamagnetism_terms[{{2,-1}}]) > tolerance) { hamiltonianmatrix += interaction_diamagnetism[{{2,-1}}]*diamagnetism_terms[{{2,-1}}]*std::sqrt(3);
}
    if (diamagnetism && std::abs(diamagnetism_terms[{{2,+2}}]) > tolerance) { hamiltonianmatrix -= interaction_diamagnetism[{{2,+2}}]*diamagnetism_terms[{{2,+2}}]*std::sqrt(1.5);
}
    if (diamagnetism && std::abs(diamagnetism_terms[{{2,-2}}]) > tolerance) { hamiltonianmatrix -= interaction_diamagnetism[{{2,-2}}]*diamagnetism_terms[{{2,-2}}]*std::sqrt(1.5);
}
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to transform the interaction /////
////////////////////////////////////////////////////////////////////

void SystemOne::transformInteraction(const eigen_sparse_t &transformator)  {
    for (auto &entry : interaction_efield) { entry.second = transformator.adjoint()*entry.second*transformator;
}
    for (auto &entry : interaction_bfield) { entry.second = transformator.adjoint()*entry.second*transformator;
}
    for (auto &entry : interaction_diamagnetism) { entry.second = transformator.adjoint()*entry.second*transformator;
}
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to delete the interaction ////////
////////////////////////////////////////////////////////////////////

void SystemOne::deleteInteraction()  {
    interaction_efield.clear();
    interaction_bfield.clear();
    interaction_diamagnetism.clear();
}


////////////////////////////////////////////////////////////////////
/// Methods that allows base class to rotate states ////////////////
////////////////////////////////////////////////////////////////////

eigen_sparse_t SystemOne::rotateStates(const std::vector<size_t> &states_indices, double alpha, double beta, double gamma) {
    // Initialize Wigner D matrix
    WignerD wigner;

    // Rotate state
    std::vector<eigen_triplet_t> states_rotated_triplets;
    states_rotated_triplets.reserve(std::min(static_cast<size_t>(10), states.size()) * states_indices.size()); // TODO std::min( 2*jmax+1, states.size() ) * states_indices.size()

    size_t current = 0;
    for (auto const &idx: states_indices) {
        this->addRotated(states[idx].state, current++, states_rotated_triplets, wigner, alpha, beta, gamma);
    }

    eigen_sparse_t states_rotated(states.size(), states_indices.size());
    states_rotated.setFromTriplets(states_rotated_triplets.begin(), states_rotated_triplets.end());
    states_rotated_triplets.clear();

    return states_rotated;
}

eigen_sparse_t SystemOne::buildStaterotator(double alpha, double beta, double gamma) {
    // Initialize Wigner D matrix
    WignerD wigner;

    // Build rotator
    std::vector<eigen_triplet_t> rotator_triplets;
    rotator_triplets.reserve(std::min(static_cast<size_t>(10), states.size()) * states.size()); // TODO std::min( 2*jmax+1, states.size() ) * states.size()

    for (auto const &entry: states) {
        this->addRotated(entry.state, entry.idx, rotator_triplets, wigner, alpha, beta, gamma);
    }

    eigen_sparse_t rotator(states.size(), states.size());
    rotator.setFromTriplets(rotator_triplets.begin(), rotator_triplets.end());
    rotator_triplets.clear();

    return rotator;
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to combine systems ///////////////
////////////////////////////////////////////////////////////////////

void SystemOne::incorporate(SystemBase<StateOne> &system) {
    // Combine parameters
    if (species != dynamic_cast<SystemOne&>(system).species) { throw std::runtime_error("The value of the variable 'element' must be the same for both systems.");
}
    if (efield != dynamic_cast<SystemOne&>(system).efield) { throw std::runtime_error("The value of the variable 'distance' must be the same for both systems."); // implies that efield_spherical is the same, too
}
    if (bfield != dynamic_cast<SystemOne&>(system).bfield) { throw std::runtime_error("The value of the variable 'angle' must be the same for both systems."); // implies that bfield_spherical is the same, too
}
    if (diamagnetism != dynamic_cast<SystemOne&>(system).diamagnetism) { throw std::runtime_error("The value of the variable 'ordermax' must be the same for both systems.");
}

    // Combine symmetries
    unsigned int num_different_symmetries = 0;
    if (sym_reflection != dynamic_cast<SystemOne&>(system).sym_reflection) {
        sym_reflection = NA;
        ++num_different_symmetries;
    }
    if (!std::equal(sym_rotation.begin(), sym_rotation.end(), dynamic_cast<SystemOne&>(system).sym_rotation.begin())) {
        if (sym_rotation.count(static_cast<float>(ARB))!=0 || dynamic_cast<SystemOne&>(system).sym_rotation.count(static_cast<float>(ARB))!=0) {
            sym_rotation = {static_cast<float>(ARB)};
        } else {
            sym_rotation.insert(dynamic_cast<SystemOne&>(system).sym_rotation.begin(), dynamic_cast<SystemOne&>(system).sym_rotation.end());
        }
        ++num_different_symmetries;
    }
    if (num_different_symmetries > 1) { std::cerr << "Warning: The systems differ in more than one symmetry. For the combined system, the notion of symmetries might be meaningless." << std::endl;
}

    // Clear cached interaction
    this->deleteInteraction();
}

////////////////////////////////////////////////////////////////////
/// Utility methods ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

void SystemOne::addCoefficient(const StateOne& state, const size_t &col, const scalar_t &value, std::vector<eigen_triplet_t> &coefficients_triplets) {
    auto state_iter = states.get<1>().find(state);

    size_t row;
    if (state_iter != states.get<1>().end()) {
        row = state_iter->idx;
    } else {
        row = states.size();
        states.push_back(enumerated_state<StateOne>(row, state));
    }

    coefficients_triplets.emplace_back(row, col, value);
}

void SystemOne::changeToSphericalbasis(std::array<double, 3> field, std::unordered_map<int, double> &field_spherical) {
    if(field[1] != 0) {
        throw std::runtime_error("For fields with non-zero y-coordinates, a complex data type is needed.");
    }
    field_spherical[1] = -field[0]/std::sqrt(2);
    field_spherical[-1] = field[0]/std::sqrt(2);
    field_spherical[0] = field[2];
}

void SystemOne::changeToSphericalbasis(std::array<double, 3> field, std::unordered_map<int, std::complex<double> > &field_spherical) {
    field_spherical[1] = std::complex<double>(-field[0]/std::sqrt(2),-field[1]/std::sqrt(2));
    field_spherical[-1] = std::complex<double>(field[0]/std::sqrt(2),-field[1]/std::sqrt(2));
    field_spherical[0] = std::complex<double>(field[2],0);
}

void SystemOne::addTriplet(std::vector<eigen_triplet_t> &triplets, const size_t r_idx, const size_t c_idx, const scalar_t val) {
    triplets.emplace_back(r_idx, c_idx, val);
    //if (r_idx != c_idx) triplets.push_back(eigen_triplet_t(c_idx, r_idx, this->conjugate(val))); // triangular matrix is not sufficient because of basis change // TODO with interaction_bfield one sometimes has to add a minus sign instead of taking the conjugate
}

void SystemOne::rotateVector(std::array<double, 3> &field, std::array<double, 3> &to_z_axis, std::array<double, 3> &to_y_axis) {
    auto field_mapped = Eigen::Map<Eigen::Matrix<double,3,1>>(&field[0]);

    if (field_mapped.norm() != 0) {
        Eigen::Matrix<double,3,3> rotator = this->buildRotator(to_z_axis, to_y_axis);
        field_mapped = rotator.transpose()*field_mapped;
    }
}

void SystemOne::rotateVector(std::array<double, 3> &field, double alpha, double beta, double gamma) {
    auto field_mapped = Eigen::Map<Eigen::Matrix<double,3,1>>(&field[0]);

    if (field_mapped.norm() != 0) {
        Eigen::Matrix<double,3,3> rotator = this->buildRotator(alpha, beta, gamma);
        field_mapped = rotator.transpose()*field_mapped;
    }
}

bool SystemOne::isRefelectionAndRotationCompatible() {
    if (sym_rotation.count(static_cast<float>(ARB))!=0 || sym_reflection == NA) { return true;
}

    for (const auto& s: sym_rotation) {
       if (sym_rotation.count(-s)==0) { return false;
}
    }

    return true;
}

