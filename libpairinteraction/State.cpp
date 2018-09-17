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

#include "State.h"
#include "QuantumDefect.h"
#include "dtypes.h"

#include <array>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>

namespace {

inline boost::variant<char, int> getMomentumLabel(int l) {
    switch (l) {
    case 0:
        return 'S';
    case 1:
        return 'P';
    case 2:
        return 'D';
    case 3:
        return 'F';
    case 4:
        return 'G';
    case 5:
        return 'H';
    case 6:
        return 'I';
    default:
        return l;
    }
}

} // namespace

std::ostream &operator<<(std::ostream &out, const State &state) {
    state.printState(out);
    return out;
} // TODO remove

////////////////////////////////////////////////////////////////////
/// Implementation of internally used state classes ////////////////
////////////////////////////////////////////////////////////////////

// StateInternalBase
std::ostream &operator<<(std::ostream &out, const StateInternalBase &state) {
    state.printState(out);
    return out;
}

bool StateInternalBase::operator==(StateInternalBase const &rhs) const {
    return (typeid(*this) == typeid(rhs)) && operatorEqual(rhs);
}

bool StateInternalBase::operator^(StateInternalBase const &rhs) const { // subset
    return (typeid(*this) == typeid(rhs)) && operatorSubset(rhs);
}

bool StateInternalBase::operator!=(StateInternalBase const &rhs) const {
    return (typeid(*this) != typeid(rhs)) || operatorUnequal(rhs);
}

bool StateInternalBase::operator<(const StateInternalBase &rhs) const {
    if (typeid(*this) != typeid(rhs)) {
        return std::type_index(typeid(*this)) < std::type_index(typeid(rhs));
    } else {
        return operatorLess(rhs);
    }
}

// StateInternalFinestructure
StateInternalFinestructure::StateInternalFinestructure(std::string species, int n, int l, float j,
                                                       float m)
    : species(std::move(species)), n(n), l(l), j(j), m(m) {
    this->analyzeSpecies();
}

bool StateInternalFinestructure::isArtificial() const { return false; }
const std::string &StateInternalFinestructure::getSpecies() const { return species; }
const std::string &StateInternalFinestructure::getElement() const { return element; }
const int &StateInternalFinestructure::getN() const { return n; }
const int &StateInternalFinestructure::getL() const { return l; }
const float &StateInternalFinestructure::getJ() const { return j; }
const float &StateInternalFinestructure::getM() const { return m; }
const float &StateInternalFinestructure::getS() const { return s; }
double StateInternalFinestructure::getEnergy() const { return energy_level(species, n, l, j); }
double StateInternalFinestructure::getNStar() const { return nstar(species, n, l, j); }

void StateInternalFinestructure::analyzeSpecies() {
    s = 0.5;
    element = species;
    if (std::isdigit(species.back()) != 0) {
        s = ((species.back() - '0') - 1) / 2.;
        element = species.substr(0, species.size() - 1);
    }
}

void StateInternalFinestructure::printState(std::ostream &out) const {
    out << species << ", " << n << " " << getMomentumLabel(l) << "_";
    if (std::ceil(j) == j) {
        out << j << ", ";
        out << "mj=" << m;
    } else {
        out << 2 * j << "/2, ";
        out << "mj=" << 2 * m << "/2";
    }
}

bool StateInternalFinestructure::operatorEqual(StateInternalBase const &rhsb) const {
    auto rhs = dynamic_cast<const StateInternalFinestructure &>(rhsb);
    return (species == rhs.species) && (n == rhs.n) && (l == rhs.l) && (j == rhs.j) && (m == rhs.m);
}
bool StateInternalFinestructure::operatorSubset(StateInternalBase const &rhsb) const { // subset
    auto rhs = dynamic_cast<const StateInternalFinestructure &>(rhsb);
    return (species == rhs.species) && (rhs.n == ARB || n == rhs.n) &&
        (rhs.l == ARB || l == rhs.l) && (rhs.j == ARB || j == rhs.j) &&
        (rhs.m == ARB || m == rhs.m);
}
bool StateInternalFinestructure::operatorUnequal(StateInternalBase const &rhsb) const {
    auto rhs = dynamic_cast<const StateInternalFinestructure &>(rhsb);
    return (species != rhs.species) || (n != rhs.n) || (l != rhs.l) || (j != rhs.j) || (m != rhs.m);
}
bool StateInternalFinestructure::operatorLess(const StateInternalBase &rhsb) const {
    auto rhs = dynamic_cast<const StateInternalFinestructure &>(rhsb);
    return (
        (species < rhs.species) ||
        ((species == rhs.species) &&
         ((n < rhs.n) ||
          ((n == rhs.n) &&
           ((l < rhs.l) || ((l == rhs.l) && ((j < rhs.j) || ((j == rhs.j) && (m < rhs.m)))))))));
}

// StateInternalArtificial
StateInternalArtificial::StateInternalArtificial(std::string label) : label(std::move(label)) {}

bool StateInternalArtificial::isArtificial() const { return true; }
const std::string &StateInternalArtificial::getLabel() const { return label; }

void StateInternalArtificial::printState(std::ostream &out) const { out << label; }

bool StateInternalArtificial::operatorEqual(StateInternalBase const &rhsb) const {
    auto rhs = dynamic_cast<const StateInternalArtificial &>(rhsb);
    return label == rhs.label;
}
bool StateInternalArtificial::operatorSubset(StateInternalBase const &rhsb) const { // subset
    auto rhs = dynamic_cast<const StateInternalArtificial &>(rhsb);
    return label == rhs.label;
}
bool StateInternalArtificial::operatorUnequal(StateInternalBase const &rhsb) const {
    auto rhs = dynamic_cast<const StateInternalArtificial &>(rhsb);
    return label != rhs.label;
}
bool StateInternalArtificial::operatorLess(const StateInternalBase &rhsb) const {
    auto rhs = dynamic_cast<const StateInternalArtificial &>(rhsb);
    return label < rhs.label;
}

////////////////////////////////////////////////////////////////////
/// Implementation of one-atom state ///////////////////////////////
////////////////////////////////////////////////////////////////////

StateOne::StateOne(std::string species, int n, int l, float j, float m)
    : state_ptr(std::make_shared<StateInternalFinestructure>(species, n, l, j, m)) {}

StateOne::StateOne(std::string label)
    : state_ptr(std::make_shared<StateInternalArtificial>(label)) {}

// Method for printing the state
std::ostream &operator<<(std::ostream &out, const StateOne &state) {
    out << "|" << *state.state_ptr << ">";
    return out;
}

// Getters
const int &StateOne::getN() const {
    return this->castState<StateInternalFinestructure>(state_ptr)->getN();
}
const int &StateOne::getL() const {
    return this->castState<StateInternalFinestructure>(state_ptr)->getL();
}
const float &StateOne::getJ() const {
    return this->castState<StateInternalFinestructure>(state_ptr)->getJ();
}
const float &StateOne::getM() const {
    return this->castState<StateInternalFinestructure>(state_ptr)->getM();
}
const float &StateOne::getS() const {
    return this->castState<StateInternalFinestructure>(state_ptr)->getS();
}
const std::string &StateOne::getSpecies() const {
    return this->castState<StateInternalFinestructure>(state_ptr)->getSpecies();
}
const std::string &StateOne::getElement() const {
    return this->castState<StateInternalFinestructure>(state_ptr)->getElement();
}
double StateOne::getEnergy() const {
    return this->castState<StateInternalFinestructure>(state_ptr)->getEnergy();
}
double StateOne::getNStar() const {
    return this->castState<StateInternalFinestructure>(state_ptr)->getNStar();
}
const std::string &StateOne::getLabel() const {
    return this->castState<StateInternalArtificial>(state_ptr)->getLabel();
}

// Comparators
bool StateOne::operator==(StateOne const &rhs) const { return *state_ptr == *rhs.state_ptr; }
bool StateOne::operator^(StateOne const &rhs) const { return *state_ptr ^ *rhs.state_ptr; }
bool StateOne::operator!=(StateOne const &rhs) const { return *state_ptr != *rhs.state_ptr; }
bool StateOne::operator<(const StateOne &rhs) const { return *state_ptr < *rhs.state_ptr; }

////////////////////////////////////////////////////////////////////
/// Implementation of StateTwo /////////////////////////////////////
////////////////////////////////////////////////////////////////////

// Constructors
StateTwo::StateTwo() : species({{"", ""}}), n({{0, 0}}), l({{0, 0}}), j({{0, 0}}), m({{0, 0}}) {
    this->analyzeSpecies();
}
StateTwo::StateTwo(std::array<std::string, 2> species, std::array<int, 2> n, std::array<int, 2> l,
                   std::array<float, 2> j, std::array<float, 2> m)
    : species(std::move(species)), n(n), l(l), j(j), m(m) {
    this->analyzeSpecies();
}
StateTwo::StateTwo(const StateOne &s1, const StateOne &s2)
    : species({{s1.getSpecies(), s2.getSpecies()}}), n({{s1.getN(), s2.getN()}}),
      l({{s1.getL(), s2.getL()}}), j({{s1.getJ(), s2.getJ()}}), m({{s1.getM(), s2.getM()}}) {
    this->analyzeSpecies();
}

// Method for printing the state
void StateTwo::printState(std::ostream &out) const {
    out << "|";
    for (size_t i = 0; i < 2; ++i) {
        out << species[i] << ", ";

        out << n[i] << " ";

        out << getMomentumLabel(l[i]) << "_";

        if (std::ceil(j[i]) == j[i]) {
            out << j[i] << ", ";
            out << "mj=" << m[i] << ">";
        } else {
            out << 2 * j[i] << "/2, ";
            out << "mj=" << 2 * m[i] << "/2";
        }

        if (i == 0) {
            out << "; ";
        }
    }
    out << ">";
}

// Comparators
bool StateTwo::operator==(const StateTwo &rhs) const {
    return (n[0] == rhs.n[0]) && (l[0] == rhs.l[0]) && (j[0] == rhs.j[0]) && (m[0] == rhs.m[0]) &&
        (n[1] == rhs.n[1]) && (l[1] == rhs.l[1]) && (j[1] == rhs.j[1]) && (m[1] == rhs.m[1]);
}

bool StateTwo::
operator^(const StateTwo &rhs) const { // subset // TODO is there a better operator to use?
    return (rhs.n[0] == ARB || n[0] == rhs.n[0]) && (rhs.l[0] == ARB || l[0] == rhs.l[0]) &&
        (rhs.j[0] == ARB || j[0] == rhs.j[0]) && (rhs.m[0] == ARB || m[0] == rhs.m[0]) &&
        (rhs.n[1] == ARB || n[1] == rhs.n[1]) && (rhs.l[1] == ARB || l[1] == rhs.l[1]) &&
        (rhs.j[1] == ARB || j[1] == rhs.j[1]) && (rhs.m[1] == ARB || m[1] == rhs.m[1]);
}

bool StateTwo::operator!=(const StateTwo &rhs) const {
    return (n[0] != rhs.n[0]) || (l[0] != rhs.l[0]) || (j[0] != rhs.j[0]) || (m[0] != rhs.m[0]) ||
        (n[1] != rhs.n[1]) || (l[1] != rhs.l[1]) || (j[1] != rhs.j[1]) || (m[1] != rhs.m[1]);
}

bool StateTwo::operator<(const StateTwo &rhs) const {
    return ((this->getFirstState() < rhs.getFirstState()) ||
            ((this->getFirstState() == rhs.getFirstState()) &&
             (this->getSecondState() < rhs.getSecondState())));
}

// Getters
double StateTwo::getEnergy() const {
    return this->getFirstState().getEnergy() + this->getSecondState().getEnergy();
}
std::array<double, 2> StateTwo::getNStar() const {
    return {{this->getFirstState().getNStar(), this->getSecondState().getNStar()}};
}
std::array<std::string, 2> StateTwo::getSpecies() const { return species; }
std::array<std::string, 2> StateTwo::getElement() const { return element; }
std::array<int, 2> StateTwo::getN() const { return n; }
std::array<int, 2> StateTwo::getL() const { return l; }
std::array<float, 2> StateTwo::getJ() const { return j; }
std::array<float, 2> StateTwo::getM() const { return m; }
std::array<float, 2> StateTwo::getS() const { return s; }
StateOne StateTwo::getFirstState() const { return StateOne(species[0], n[0], l[0], j[0], m[0]); }
StateOne StateTwo::getSecondState() const { return StateOne(species[1], n[1], l[1], j[1], m[1]); }

// Utility methods
void StateTwo::analyzeSpecies() {
    for (size_t i = 0; i < 2; ++i) {
        s[i] = 0.5;
        element[i] = species[i];

        if (std::isdigit(species[i].back()) != 0) {
            s[i] = ((species[i].back() - '0') - 1) / 2.;
            element[i] = species[i].substr(0, species[i].size() - 1);
        }
    }
}
