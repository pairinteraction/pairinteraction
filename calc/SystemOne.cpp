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

SystemOne::SystemOne(std::string const& element, std::string cachedir)
    : SystemBase(cachedir), efield({{0,0,0}}), bfield({{0,0,0}}), diamagnetism(false), element(element)
{}

const std::string& SystemOne::getElement() const {
    return element;
}

void SystemOne::setEfield(std::array<double, 3> field) {
    efield = field;
}

void SystemOne::setBfield(std::array<double, 3> field) {
    bfield = field;
}

void SystemOne::setDiamagnetism(bool enable) {
    diamagnetism = enable;
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to initialize Basis //////////////
////////////////////////////////////////////////////////////////////

void SystemOne::initializeBasis()
{
    // TODO check whether specified basis is finite

    ////////////////////////////////////////////////////////////////////
    /// Build one atom states //////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    // TODO consider symmetries

    size_t idx = 0;
    std::vector<eigen_triplet_t> coefficients_triplets; // TODO reserve states, coefficients_triplets, hamiltonianmatrix_triplets
    std::vector<eigen_triplet_t> hamiltonianmatrix_triplets;
    std::set<int> range_adapted_n, range_adapted_l;
    std::set<float> range_adapted_j, range_adapted_m;

    if (range_n.empty()) {
        range_adapted_n = std::set<int>({}); // TODO if empty, calculate the range via the energies
    } else {
        range_adapted_n = range_n;
    }
    for (auto n : range_adapted_n) {

        if (range_l.empty()) {
            range_adapted_l.clear();
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
                if (!checkIsEnergyValid(energy)) continue;

                if (range_m.empty()) {
                    range_adapted_m.clear();
                    for (float m = -j; m < j; ++m) {
                        range_adapted_m.insert(m);
                    }
                } else {
                    range_adapted_m = range_m;
                }
                for (auto m : range_adapted_m) {
                    if (std::fabs(m) > j) continue;

                    states.push_back(StateOne(element,n,l,j,m));
                    hamiltonianmatrix_triplets.push_back(eigen_triplet_t(idx,idx,energy));
                    coefficients_triplets.push_back(eigen_triplet_t(idx,idx,1)); // TODO take into account symmetries

                    ++idx;
                }
            }
        }
    }

    // Build data
    states.shrink_to_fit();

    coefficients.resize(idx,idx);
    coefficients.setFromTriplets(coefficients_triplets.begin(), coefficients_triplets.end());
    coefficients_triplets.clear();

    hamiltonianmatrix.resize(idx,idx);
    hamiltonianmatrix.setFromTriplets(hamiltonianmatrix_triplets.begin(), hamiltonianmatrix_triplets.end());
    hamiltonianmatrix_triplets.clear();
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to initialize Hamiltonian helpers
////////////////////////////////////////////////////////////////////

void SystemOne::initializeHamiltonianhelpers() {
    ////////////////////////////////////////////////////////////////////
    /// Prepare the calculation of the Hamiltonian helpers /////////////
    ////////////////////////////////////////////////////////////////////

    // Transform the electromagnetic fields into spherical coordinates
    std::unordered_map<int, scalar_t>  efield_spherical, bfield_spherical;
    this->changeToSphericalbasis(efield, efield_spherical);
    this->changeToSphericalbasis(efield, bfield_spherical);

    // Check if something to do and delete Hamiltonian helpers if they are not needed
    std::vector<int> erange, brange;
    for (auto entry : efield_spherical) {
        if (std::abs(entry.second) == 0) {
            hamiltonianhelper_efield.erase(-entry.first);
        } else if (hamiltonianhelper_efield.find(-entry.first) == hamiltonianhelper_efield.end()) {
            erange.push_back(-entry.first);
        }
    }
    for (auto entry : bfield_spherical) {
        if (std::abs(entry.second) == 0) {
            hamiltonianhelper_bfield.erase(-entry.first);
        } else if (hamiltonianhelper_bfield.find(-entry.first) == hamiltonianhelper_bfield.end()) {
            brange.push_back(-entry.first);
        }
    }

    if (erange.empty() && brange.empty()) return;

    // TODO add operators for diamagnetism !!!

    // Precalculate matrix elements
    MatrixElements matrix_elements(element, (cachedir / "cache_elements.db").string());
    for (int i : erange) matrix_elements.precalculateElectricMomentum(states, i);
    for (int i : brange) matrix_elements.precalculateMagneticMomentum(states, i);

    // TODO add operators for diamagnetism !!!

    ////////////////////////////////////////////////////////////////////
    /// Generate the Hamiltonian helpers in the canonical basis ////////
    ////////////////////////////////////////////////////////////////////

    std::unordered_map<int, std::vector<eigen_triplet_t>> hamiltonianhelper_efield_triplets; // TODO reserve
    std::unordered_map<int, std::vector<eigen_triplet_t>> hamiltonianhelper_bfield_triplets;// TODO reserve
    std::unordered_map<std::array<int, 2>, std::vector<eigen_triplet_t>> hamiltonianhelper_diamagnetism_triplets; // TODO reserve

    for (size_t col=0; col<states.size(); ++col) { // TODO parallelization
        for (size_t row=0; row<states.size(); ++row) {
            //if (row < col) continue; // TODO use this restriction and construct the total matrix (that is needed in order to be transformable) from the diagonal matrix afterwards

            const StateOne &state_row = states[row];
            const StateOne &state_col = states[col];

            for (int i : erange) {
                if (selectionRulesMultipole(state_row, state_col, 1, i)) {
                    hamiltonianhelper_efield_triplets[i].push_back(eigen_triplet_t(row, col, matrix_elements.getElectricMomentum(state_row, state_col)));
                    break;
                }
            }

            for (int i : brange) {
                if (selectionRulesMomentum(state_row, state_col, i)) {
                    hamiltonianhelper_bfield_triplets[i].push_back(eigen_triplet_t(row, col, matrix_elements.getMagneticMomentum(state_row, state_col)));
                    break;
                }
            }

            // TODO add operators for diamagnetism !!!
        }
    }

    for (auto &triplets : hamiltonianhelper_efield_triplets) {
        hamiltonianhelper_efield[triplets.first].resize(states.size(),states.size());
        hamiltonianhelper_efield[triplets.first].setFromTriplets(triplets.second.begin(), triplets.second.end());
        triplets.second.clear();
    }

    for (auto &triplets : hamiltonianhelper_bfield_triplets) {
        hamiltonianhelper_bfield[triplets.first].resize(states.size(),states.size());
        hamiltonianhelper_bfield[triplets.first].setFromTriplets(triplets.second.begin(), triplets.second.end());
        triplets.second.clear();
    }

    for (auto &triplets : hamiltonianhelper_diamagnetism_triplets) {
        hamiltonianhelper_diamagnetism[triplets.first].resize(states.size(),states.size());
        hamiltonianhelper_diamagnetism[triplets.first].setFromTriplets(triplets.second.begin(), triplets.second.end());
        triplets.second.clear();
    }

    ////////////////////////////////////////////////////////////////////
    /// Transform the Hamiltonian helpers to the used basis ////////////
    ////////////////////////////////////////////////////////////////////

    for (auto &helper : hamiltonianhelper_efield) helper.second = coefficients.adjoint()*helper.second*coefficients;
    for (auto &helper : hamiltonianhelper_bfield) helper.second = coefficients.adjoint()*helper.second*coefficients;
    for (auto &helper : hamiltonianhelper_diamagnetism) helper.second = coefficients.adjoint()*helper.second*coefficients;
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to construct Hamiltonian /////////
////////////////////////////////////////////////////////////////////

void SystemOne::initializeHamiltonian() {
    // Update the Hamiltonian helpers // TODO initializeHamiltonianhelpers needs not to be called by SystemBase explicitely
    this->initializeHamiltonianhelpers(); // TODO rename to buildHamiltonianhelpers()

    // Transform the electromagnetic fields into spherical coordinates
    std::unordered_map<int, scalar_t> efield_spherical, bfield_spherical;
    this->changeToSphericalbasis(efield, efield_spherical);
    this->changeToSphericalbasis(efield, bfield_spherical);

    // Build the total Hamiltonian
    if (std::abs(efield_spherical[0]) != 0) hamiltonianmatrix -= hamiltonianhelper_efield[0]*efield_spherical[0];
    if (std::abs(efield_spherical[-1]) != 0) hamiltonianmatrix += hamiltonianhelper_efield[1]*efield_spherical[-1];
    if (std::abs(efield_spherical[1]) != 0) hamiltonianmatrix += hamiltonianhelper_efield[-1]*efield_spherical[1];
    if (std::abs(bfield_spherical[0]) != 0) hamiltonianmatrix += hamiltonianhelper_bfield[0]*bfield_spherical[0];
    if (std::abs(bfield_spherical[-1]) != 0) hamiltonianmatrix -= hamiltonianhelper_bfield[1]*bfield_spherical[-1];
    if (std::abs(bfield_spherical[1]) != 0) hamiltonianmatrix -= hamiltonianhelper_bfield[-1]*bfield_spherical[1];
    // TODO add operators for diamagnetism !!!

    // Delete the Hamiltonian helpers as they are no longer needed
    hamiltonianhelper_efield.clear();
    hamiltonianhelper_bfield.clear();
    hamiltonianhelper_diamagnetism.clear();

    // TODO throw an error if efield etc. is changed after initializeHamiltonian was called, or prevent it from beeing changed by introducing a SystemGenerator
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to transform Hamiltonian helpers /
////////////////////////////////////////////////////////////////////

void SystemOne::transformHamiltonianhelpers(const eigen_sparse_t &transformator)  {
    for (auto &helper : hamiltonianhelper_efield) helper.second = transformator.adjoint()*helper.second*transformator;
    for (auto &helper : hamiltonianhelper_bfield) helper.second = transformator.adjoint()*helper.second*transformator;
    for (auto &helper : hamiltonianhelper_diamagnetism) helper.second = transformator.adjoint()*helper.second*transformator;
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

