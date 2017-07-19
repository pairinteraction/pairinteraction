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
#include "MatrixElements.h"

#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>
#include <unordered_set>
#include <type_traits>

SystemOne::SystemOne(std::wstring const& element, std::wstring cachedir)
    : SystemBase(cachedir), efield({{0,0,0}}), bfield({{0,0,0}}), diamagnetism(false), element(element){
}

SystemOne::SystemOne(std::wstring const& element, std::wstring cachedir, bool memory_saving)
    : SystemBase(cachedir, memory_saving), efield({{0,0,0}}), bfield({{0,0,0}}), diamagnetism(false), element(element){
}

SystemOne::SystemOne(std::wstring const& element)
    : SystemBase(), efield({{0,0,0}}), bfield({{0,0,0}}), diamagnetism(false), element(element){
}

SystemOne::SystemOne(std::wstring const& element, bool memory_saving)
    : SystemBase(memory_saving), efield({{0,0,0}}), bfield({{0,0,0}}), diamagnetism(false), element(element){
}

const std::wstring& SystemOne::getElement() const {
    return element;
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

    // TODO consider symmetries and check whether they are applicable

    size_t idx = 0;
    std::vector<eigen_triplet_t> coefficients_triplets; // TODO reserve states, coefficients_triplets, hamiltonianmatrix_triplets
    std::vector<eigen_triplet_t> hamiltonianmatrix_triplets;
    std::set<int> range_adapted_n, range_adapted_l;
    std::set<float> range_adapted_j, range_adapted_m;

    if (range_n.empty()) {
        throw std::runtime_error( "The calculation of range_n via energy restrictions is not yet implemented." ); // TODO
    } else {
        range_adapted_n = range_n;
    }
    for (auto n : range_adapted_n) {

        if (range_l.empty()) {
            this->range(range_adapted_l, 0, n-1);
        } else {
            range_adapted_l = range_l;
        }
        for (auto l : range_adapted_l) {
            if (l > n-1) continue;

            if (range_j.empty()) {
                range_adapted_j = {std::fabs(l-0.5f), l+0.5f};
            } else {
                range_adapted_j = range_j;
            }
            for (auto j : range_adapted_j) {
                if (std::fabs(j-l) != 0.5) continue;

                double energy = StateOne(element,n,l,j,0.5).getEnergy();
                if (!checkIsEnergyValid(energy)) continue;

                if (range_m.empty()) {
                    this->range(range_adapted_m, -j, j);
                } else {
                    range_adapted_m = range_m;
                }
                for (auto m : range_adapted_m) {
                    if (std::fabs(m) > j) continue;

                    states.push_back(enumerated_state<StateOne>(idx, StateOne(element,n,l,j,m)));
                    hamiltonianmatrix_triplets.push_back(eigen_triplet_t(idx,idx,energy));
                    coefficients_triplets.push_back(eigen_triplet_t(idx,idx,1));

                    ++idx;
                }
            }
        }
    }

    // Build data
    coefficients.resize(idx,idx);
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

    if (erange.empty() && brange.empty() && drange.empty()) return;

    // Precalculate matrix elements
    std::string matrixelementsdir = "";
    if (!cachedir.empty()) matrixelementsdir = (cachedir / "cache_elements.db").string(); // TODO do this in the MatrixElements class, just pass cachedir as an argument to the constructor

    std::string tmp(element.begin(), element.end()); // TODO think of a better solution
    MatrixElements matrixelements(tmp, matrixelementsdir);
    auto states_converted = this->getStates(); // TODO remove
    for (const auto &i : erange) matrixelements.precalculateElectricMomentum(states_converted, i);
    for (const auto &i : brange) matrixelements.precalculateMagneticMomentum(states_converted, i);
    for (const auto &i : drange) matrixelements.precalculateDiamagnetism(states_converted, i[0], i[1]);

    ////////////////////////////////////////////////////////////////////
    /// Calculate the interaction in the canonical basis ///////////////
    ////////////////////////////////////////////////////////////////////

    std::unordered_map<int, std::vector<eigen_triplet_t>> interaction_efield_triplets; // TODO reserve
    std::unordered_map<int, std::vector<eigen_triplet_t>> interaction_bfield_triplets;// TODO reserve
    std::unordered_map<std::array<int, 2>, std::vector<eigen_triplet_t>> interaction_diamagnetism_triplets; // TODO reserve

    // Loop over column entries
    for (const auto &c: states) { // TODO parallelization

        if (c.state.element.empty()) continue; // TODO artifical states TODO [dummystates]

        // Loop over row entries
        for (const auto &r: states) {

            if (r.state.element.empty()) continue; // TODO artifical states TODO [dummystates]
            if (r.idx < c.idx) continue;

            // E-field interaction
            for (const auto &i : erange) {
                if (selectionRulesMultipole(r.state, c.state, 1, i)) {
                    scalar_t value = matrixelements.getElectricMomentum(r.state, c.state);
                    this->addTriplet(interaction_efield_triplets[i], r.idx, c.idx, value);
                    break; // because for the other operators, the selection rule for the magnetic quantum numbers will not be fulfilled
                }
            }

            // B-field interaction
            for (const auto &i : brange) {
                if (selectionRulesMomentum(r.state, c.state, i)) {
                    scalar_t value = matrixelements.getMagneticMomentum(r.state, c.state);
                    this->addTriplet(interaction_bfield_triplets[i], r.idx, c.idx, value);
                    break; // because for the other operators, the selection rule for the magnetic quantum numbers will not be fulfilled
                }
            }

            // Diamagnetic interaction
            for (const auto &i : drange) {
                if (selectionRulesMultipole(r.state, c.state, i[0], i[1])) {
                    scalar_t value = matrixelements.getDiamagnetism(r.state, c.state, i[0]);
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

    if (std::abs(efield_spherical[+0]) > tolerance) hamiltonianmatrix -= interaction_efield[+0]*efield_spherical[+0];
    if (std::abs(efield_spherical[-1]) > tolerance) hamiltonianmatrix += interaction_efield[+1]*efield_spherical[-1];
    if (std::abs(efield_spherical[+1]) > tolerance) hamiltonianmatrix += interaction_efield[-1]*efield_spherical[+1];
    if (std::abs(bfield_spherical[+0]) > tolerance) hamiltonianmatrix += interaction_bfield[+0]*bfield_spherical[+0];
    if (std::abs(bfield_spherical[-1]) > tolerance) hamiltonianmatrix -= interaction_bfield[+1]*bfield_spherical[-1];
    if (std::abs(bfield_spherical[+1]) > tolerance) hamiltonianmatrix -= interaction_bfield[-1]*bfield_spherical[+1];

    if (diamagnetism && std::abs(diamagnetism_terms[{{0,+0}}]) > tolerance) hamiltonianmatrix += interaction_diamagnetism[{{0,+0}}]*diamagnetism_terms[{{0,+0}}];
    if (diamagnetism && std::abs(diamagnetism_terms[{{2,+0}}]) > tolerance) hamiltonianmatrix -= interaction_diamagnetism[{{2,+0}}]*diamagnetism_terms[{{2,+0}}];
    if (diamagnetism && std::abs(diamagnetism_terms[{{2,+1}}]) > tolerance) hamiltonianmatrix += interaction_diamagnetism[{{2,+1}}]*diamagnetism_terms[{{2,+1}}]*std::sqrt(3);
    if (diamagnetism && std::abs(diamagnetism_terms[{{2,-1}}]) > tolerance) hamiltonianmatrix += interaction_diamagnetism[{{2,-1}}]*diamagnetism_terms[{{2,-1}}]*std::sqrt(3);
    if (diamagnetism && std::abs(diamagnetism_terms[{{2,+2}}]) > tolerance) hamiltonianmatrix -= interaction_diamagnetism[{{2,+2}}]*diamagnetism_terms[{{2,+2}}]*std::sqrt(1.5);
    if (diamagnetism && std::abs(diamagnetism_terms[{{2,-2}}]) > tolerance) hamiltonianmatrix -= interaction_diamagnetism[{{2,-2}}]*diamagnetism_terms[{{2,-2}}]*std::sqrt(1.5);
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to transform the interaction /////
////////////////////////////////////////////////////////////////////

void SystemOne::transformInteraction(const eigen_sparse_t &transformator)  {
    for (auto &entry : interaction_efield) entry.second = transformator.adjoint()*entry.second*transformator;
    for (auto &entry : interaction_bfield) entry.second = transformator.adjoint()*entry.second*transformator;
    for (auto &entry : interaction_diamagnetism) entry.second = transformator.adjoint()*entry.second*transformator;
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

eigen_sparse_t SystemOne::rotateState(const StateOne &state, double alpha, double beta, double gamma) {
    // Initialize Wigner D matrix
    WignerD wigner;

    // Rotate state
    std::vector<eigen_triplet_t> state_rotated_triplets;
    state_rotated_triplets.reserve(std::min(static_cast<size_t>(2*state.j+1), states.size()));

    this->addRotated(state, 0, state_rotated_triplets, wigner, alpha, beta, gamma);

    eigen_sparse_t state_rotated(states.size(),1);
    state_rotated.setFromTriplets(state_rotated_triplets.begin(), state_rotated_triplets.end());
    state_rotated_triplets.clear();

    return state_rotated;
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
/// Utility methods ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

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
    triplets.push_back(eigen_triplet_t(r_idx, c_idx, val));
    if (r_idx != c_idx) triplets.push_back(eigen_triplet_t(c_idx, r_idx, this->conjugate(val))); // triangular matrix is not sufficient because of basis change
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

template<>
double SystemOne::convert(const std::complex<double> &val) {
    return val.real();
}

