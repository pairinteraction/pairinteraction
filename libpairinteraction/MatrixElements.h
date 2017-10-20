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

#ifndef MATRIXELEMENTS_H
#define MATRIXELEMENTS_H

#include "dtypes.h"
#include "Basisnames.h"
#include "Wavefunction.h"
#include "wignerSymbols/include/wignerSymbols/wignerSymbols-cpp.h"

#include <string>
#include <unordered_map>
#include <memory>
#include <sstream>

bool selectionRulesMomentum(StateOne const& state1, StateOne const& state2, int q);
bool selectionRulesMomentum(StateOne const& state1, StateOne const& state2);
bool selectionRulesMultipole(StateOne const& state1, StateOne const& state2, int kappa, int q);
bool selectionRulesMultipole(StateOne const& state1, StateOne const& state2, int kappa);

class MatrixElements { // TODO use one single buffer, initialiezed at the start of the program
public:
    MatrixElements(std::string const& species, std::string const& dbname);
    MatrixElements(const Configuration& config, std::string const& species, std::string const& dbname);
    void precalculateElectricMomentum(std::shared_ptr<const BasisnamesOne> basis_one, int q);
    void precalculateMagneticMomentum(std::shared_ptr<const BasisnamesOne> basis_one, int q);
    void precalculateDiamagnetism(std::shared_ptr<const BasisnamesOne> basis_one, int k, int q);
    void precalculateMultipole(std::shared_ptr<const BasisnamesOne> basis_one, int k);
    void precalculateRadial(std::shared_ptr<const BasisnamesOne> basis_one, int k);
    double getElectricMomentum(StateOne const& state_row, StateOne const& state_col);
    double getMagneticMomentum(StateOne const& state_row, StateOne const& state_col);
    double getDiamagnetism(StateOne const& state_row, StateOne const& state_col, int kappa);
    double getMultipole(StateOne const& state_row, StateOne const& state_col, int kappa);
    double getRadial(StateOne const& state_row, StateOne const& state_col, int kappa);

    void precalculateElectricMomentum(const std::vector<StateOne> &basis_one, int q);
    void precalculateMagneticMomentum(const std::vector<StateOne> &basis_one, int q);
    void precalculateDiamagnetism(const std::vector<StateOne> &basis_one, int k, int q);
    void precalculateMultipole(const std::vector<StateOne> &basis_one, int k);
    void precalculateRadial(const std::vector<StateOne> &basis_one, int k);
private:
    void precalculate(std::shared_ptr<const BasisnamesOne> basis_one, int kappa, int q, int kappar, bool calcMultipole, bool calcMomentum, bool calcRadial);
    double calcRadialElement(const QuantumDefect &qd1, int power, const QuantumDefect &qd2);
    std::string method;
    std::string species;
    std::string dbname;
    std::unordered_map<int,std::unordered_map<StateTwo,double>> cache_radial;
    std::unordered_map<int,std::unordered_map<StateTwo,double>> cache_angular;
    std::unordered_map<int,std::unordered_map<StateTwo,double>> cache_reduced_commutes_s;
    std::unordered_map<int,std::unordered_map<StateTwo,double>> cache_reduced_commutes_l;
    std::unordered_map<int,std::unordered_map<StateTwo,double>> cache_reduced_multipole;
    double muB; // TODO define them in constants.h
    double gS;
    double gL;

    void precalculate(const std::vector<StateOne> &basis_one, int kappa, int q, int kappar, bool calcMultipole, bool calcMomentum, bool calcRadial);
};

#endif
