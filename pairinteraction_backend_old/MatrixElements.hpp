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

#ifndef MATRIXELEMENTS_H
#define MATRIXELEMENTS_H

#include "Basisnames.hpp"
#include "StateOld.hpp"
#include "Wavefunction.hpp"

#include <wignerSymbols/wignerSymbols-cpp.h>

#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>

bool selectionRulesMomentum(StateOneOld const &state1, StateOneOld const &state2, int q);
bool selectionRulesMomentum(StateOneOld const &state1, StateOneOld const &state2);
bool selectionRulesMultipole(StateOneOld const &state1, StateOneOld const &state2, int kappa,
                             int q);
bool selectionRulesMultipole(StateOneOld const &state1, StateOneOld const &state2, int kappa);

class MatrixElements { // TODO use one single buffer, initialiezed at the start of the program
public:
    MatrixElements(std::string const &species, std::string dbname);
    MatrixElements(const Configuration &config, std::string const &species,
                   std::string const &dbname);
    void precalculateElectricMomentum(std::shared_ptr<const BasisnamesOne> const &basis_one, int q);
    void precalculateMagneticMomentum(std::shared_ptr<const BasisnamesOne> const &basis_one, int q);
    void precalculateDiamagnetism(std::shared_ptr<const BasisnamesOne> const &basis_one, int k,
                                  int q);
    void precalculateMultipole(std::shared_ptr<const BasisnamesOne> const &basis_one, int k);
    void precalculateRadial(std::shared_ptr<const BasisnamesOne> const &basis_one, int k);
    double getElectricMomentum(StateOneOld const &state_row, StateOneOld const &state_col);
    double getMagneticMomentum(StateOneOld const &state_row, StateOneOld const &state_col);
    double getDiamagnetism(StateOneOld const &state_row, StateOneOld const &state_col, int k);
    double getMultipole(StateOneOld const &state_row, StateOneOld const &state_col, int k);
    double getRadial(StateOneOld const &state_row, StateOneOld const &state_col, int k);

    void precalculateElectricMomentum(const std::vector<StateOneOld> &basis_one, int q);
    void precalculateMagneticMomentum(const std::vector<StateOneOld> &basis_one, int q);
    void precalculateDiamagnetism(const std::vector<StateOneOld> &basis_one, int k, int q);
    void precalculateMultipole(const std::vector<StateOneOld> &basis_one, int k);
    void precalculateRadial(const std::vector<StateOneOld> &basis_one, int k);

private:
    void precalculate(const std::shared_ptr<const BasisnamesOne> &basis_one, int kappa, int q,
                      int kappar, bool calcElectricMultipole, bool calcMagneticMomentum,
                      bool calcRadial);
    double calcRadialElement(const QuantumDefect &qd1, int power, const QuantumDefect &qd2);
    std::string method;
    std::string species;
    std::string dbname;
    std::unordered_map<int, std::unordered_map<StateTwoOld, double>> cache_radial;
    std::unordered_map<int, std::unordered_map<StateTwoOld, double>> cache_angular;
    std::unordered_map<int, std::unordered_map<StateTwoOld, double>> cache_reduced_commutes_s;
    std::unordered_map<int, std::unordered_map<StateTwoOld, double>> cache_reduced_commutes_l;
    std::unordered_map<int, std::unordered_map<StateTwoOld, double>> cache_reduced_multipole;
    double muB; // TODO define them in constants.h
    double gS;
    double gL;
    double s;

    void precalculate(const std::vector<StateOneOld> &basis_one, int kappa, int q, int kappar,
                      bool calcElectricMultipole, bool calcMagneticMomentum, bool calcRadial);
};

#endif
