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
}

void SystemOne::setBfield(std::array<double, 3> field) {
    this->onParameterChange();
    bfield = field;
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

                    states.push_back(StateOne(element,n,l,j,m));
                    hamiltonianmatrix_triplets.push_back(eigen_triplet_t(idx,idx,energy));
                    coefficients_triplets.push_back(eigen_triplet_t(idx,idx,1));

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
/// Method that allows base class to calculate the interaction /////
////////////////////////////////////////////////////////////////////

void SystemOne::initializeInteraction() {
    ////////////////////////////////////////////////////////////////////
    /// Prepare the calculation of the interaction /////////////////////
    ////////////////////////////////////////////////////////////////////

    // Transform the electromagnetic fields into spherical coordinates
    std::unordered_map<int, scalar_t>  efield_spherical, bfield_spherical;
    this->changeToSphericalbasis(efield, efield_spherical);
    this->changeToSphericalbasis(bfield, bfield_spherical);

    // Check if something to do
    std::vector<int> erange, brange;
    std::vector<std::array<int, 2>> drange;
    for (auto entry : efield_spherical) {
        if (std::abs(entry.second) != 0 && interaction_efield.find(-entry.first) == interaction_efield.end()) {
            erange.push_back(-entry.first);
        }
    }
    for (auto entry : bfield_spherical) {
        if (std::abs(entry.second) != 0 && interaction_bfield.find(-entry.first) == interaction_bfield.end()) {
            brange.push_back(-entry.first);
        }
    }

    if (erange.empty() && brange.empty() && drange.empty()) return;

    // TODO add operators for diamagnetism !!!

    // Precalculate matrix elements
    std::string matrixelementsdir = "";
    if (!cachedir.empty()) matrixelementsdir = (cachedir / "cache_elements.db").string(); // TODO do this in the MatrixElements class, just pass cachedir as an argument to the constructor

    MatrixElements matrixelements(std::wstring_convert<std::codecvt_utf8<wchar_t>>().to_bytes(element), matrixelementsdir);
    for (int i : erange) matrixelements.precalculateElectricMomentum(states, i);
    for (int i : brange) matrixelements.precalculateMagneticMomentum(states, i);

    // TODO add operators for diamagnetism !!!

    ////////////////////////////////////////////////////////////////////
    /// Generate the interaction in the canonical basis ////////////////
    ////////////////////////////////////////////////////////////////////

    std::unordered_map<int, std::vector<eigen_triplet_t>> interaction_efield_triplets; // TODO reserve
    std::unordered_map<int, std::vector<eigen_triplet_t>> interaction_bfield_triplets;// TODO reserve
    std::unordered_map<std::array<int, 2>, std::vector<eigen_triplet_t>> interaction_diamagnetism_triplets; // TODO reserve

    for (size_t col=0; col<states.size(); ++col) { // TODO parallelization
        const StateOne &state_col = states[col];
        if (state_col.element.empty()) continue; // TODO artifical states TODO [dummystates]

        for (size_t row=0; row<states.size(); ++row) {
            const StateOne &state_row = states[row];
            if (state_row.element.empty()) continue; // TODO artifical states TODO [dummystates]

            if (row < col) continue;

            for (int i : erange) {
                if (selectionRulesMultipole(state_row, state_col, 1, i)) {
                    scalar_t value = matrixelements.getElectricMomentum(state_row, state_col);
                    interaction_efield_triplets[i].push_back(eigen_triplet_t(row, col, value));
                    if (row != col) interaction_efield_triplets[i].push_back(eigen_triplet_t(col, row, this->conjugate(value)));
                    break;
                }
            }

            for (int i : brange) {
                if (selectionRulesMomentum(state_row, state_col, i)) {
                    scalar_t value = matrixelements.getMagneticMomentum(state_row, state_col);
                    interaction_bfield_triplets[i].push_back(eigen_triplet_t(row, col, value));
                    if (row != col) interaction_bfield_triplets[i].push_back(eigen_triplet_t(col, row, this->conjugate(value)));
                    break;
                }
            }

            // TODO add operators for diamagnetism !!!
        }
    }

    for (const auto &i : erange) {
        interaction_efield[i].resize(states.size(),states.size());
        interaction_efield[i].setFromTriplets(interaction_efield_triplets[i].begin(), interaction_efield_triplets[i].end());
        interaction_efield_triplets[i].clear();
    }

    for (const auto &i : brange) {
        interaction_bfield[i].resize(states.size(),states.size());
        interaction_bfield[i].setFromTriplets(interaction_bfield_triplets[i].begin(), interaction_bfield_triplets[i].end());
        interaction_bfield_triplets[i].clear();
    }

    for (const auto &i : drange) {
        interaction_diamagnetism[i].resize(states.size(),states.size());
        interaction_diamagnetism[i].setFromTriplets(interaction_diamagnetism_triplets[i].begin(), interaction_diamagnetism_triplets[i].end());
        interaction_diamagnetism_triplets[i].clear();
    }

    ////////////////////////////////////////////////////////////////////
    /// Transform the interaction to the used basis ////////////////////
    ////////////////////////////////////////////////////////////////////

    for (auto &entry : interaction_efield) entry.second = coefficients.adjoint()*entry.second*coefficients;
    for (auto &entry : interaction_bfield) entry.second = coefficients.adjoint()*entry.second*coefficients;
    for (auto &entry : interaction_diamagnetism) entry.second = coefficients.adjoint()*entry.second*coefficients;
}

////////////////////////////////////////////////////////////////////
/// Method that allows base class to construct Hamiltonian /////////
////////////////////////////////////////////////////////////////////

void SystemOne::addInteraction() {

    // Transform the electromagnetic fields into spherical coordinates
    std::unordered_map<int, scalar_t> efield_spherical, bfield_spherical;
    this->changeToSphericalbasis(efield, efield_spherical);
    this->changeToSphericalbasis(bfield, bfield_spherical);

    // Build the total Hamiltonian
    if (std::abs(efield_spherical[0]) != 0) hamiltonianmatrix -= interaction_efield[0]*efield_spherical[0];
    if (std::abs(efield_spherical[-1]) != 0) hamiltonianmatrix += interaction_efield[1]*efield_spherical[-1];
    if (std::abs(efield_spherical[1]) != 0) hamiltonianmatrix += interaction_efield[-1]*efield_spherical[1];
    if (std::abs(bfield_spherical[0]) != 0) hamiltonianmatrix += interaction_bfield[0]*bfield_spherical[0];
    if (std::abs(bfield_spherical[-1]) != 0) hamiltonianmatrix -= interaction_bfield[1]*bfield_spherical[-1];
    if (std::abs(bfield_spherical[1]) != 0) hamiltonianmatrix -= interaction_bfield[-1]*bfield_spherical[1];

    // TODO add operators for diamagnetism !!!
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
