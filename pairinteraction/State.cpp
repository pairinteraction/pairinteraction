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

#include "State.h"
#include "QuantumDefect.h"
#include "dtypes.h"
#include "utils.h"

#include <array>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>

namespace {

struct getMomentumLabel {
    int l;
    getMomentumLabel(int l) : l(l) {}
};

inline std::ostream &operator<<(std::ostream &os, getMomentumLabel const &l) {
    static constexpr std::array const letters = {'S', 'P', 'D', 'F', 'G', 'H', 'I'};
    if (l.l < static_cast<int>(letters.size())) {
        os << letters[l.l];
    } else {
        os << l.l;
    }
    return os;
}

} // namespace

////////////////////////////////////////////////////////////////////
/// Implementation of one-atom state ///////////////////////////////
////////////////////////////////////////////////////////////////////

StateOne::StateOne(std::string species, int n, int l, float j, float m)
    : species(std::move(species)), n(n), l(l), j(j), m(m) {
    this->analyzeSpecies();
    hashvalue = 0;
    utils::hash_combine(hashvalue, this->getSpecies());
    utils::hash_combine(hashvalue, this->getN());
    utils::hash_combine(hashvalue, this->getL());
    utils::hash_combine(hashvalue, this->getJ());
    utils::hash_combine(hashvalue, this->getM());
}

StateOne::StateOne(std::string label)
    : species(std::move(label)), element(""), n(0), l(0), j(0), m(0), s(0) {
    hashvalue = std::hash<std::string>{}(this->getLabel());
}

// Methods for printing the state
std::ostream &operator<<(std::ostream &out, const StateOne &state) {
    out << "|";
    if (state.isArtificial()) {
        out << state.getLabel();
    } else {
        out << state.getSpecies() << ", " << state.getN() << " " << getMomentumLabel(state.getL())
            << "_";
        if (std::ceil(state.getJ()) == state.getJ()) {
            out << state.getJ() << ", ";
            out << "mj=" << state.getM();
        } else {
            out << 2 * state.getJ() << "/2, ";
            out << "mj=" << 2 * state.getM() << "/2";
        }
    }
    out << ">";
    return out;
}

std::string StateOne::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

// Getters
const int &StateOne::getN() const {
    this->shouldBeArtificial(false);
    return n;
}
const int &StateOne::getL() const {
    this->shouldBeArtificial(false);
    return l;
}
const float &StateOne::getJ() const {
    this->shouldBeArtificial(false);
    return j;
}
const float &StateOne::getM() const {
    this->shouldBeArtificial(false);
    return m;
}
const float &StateOne::getS() const {
    this->shouldBeArtificial(false);
    return s;
}
const std::string &StateOne::getSpecies() const {
    this->shouldBeArtificial(false);
    return species;
}
const std::string &StateOne::getElement() const {
    this->shouldBeArtificial(false);
    return element;
}
double StateOne::getEnergy() const {
    this->shouldBeArtificial(false);
    return energy_level(species, n, l, j);
}
double StateOne::getEnergy(MatrixElementCache &cache) const {
    this->shouldBeArtificial(false);
    return energy_level(species, n, l, j, cache.getDefectDB());
}
double StateOne::getNStar() const {
    this->shouldBeArtificial(false);
    return nstar(species, n, l, j);
}
double StateOne::getNStar(MatrixElementCache &cache) const {
    this->shouldBeArtificial(false);
    return nstar(species, n, l, j, cache.getDefectDB());
}
const std::string &StateOne::getLabel() const {
    this->shouldBeArtificial(true);
    return species;
}
bool StateOne::isArtificial() const { return (n == 0); }
bool StateOne::isGeneralized() const {
    return (n == ARB) || (l == ARB) || (j == ARB) || (m == ARB);
}

const size_t &StateOne::getHash() const { return hashvalue; }

StateOne StateOne::getReflected() const {
    return StateOne(this->getSpecies(), this->getN(), this->getL(), this->getJ(), -this->getM());
}

// Comparators
bool StateOne::operator==(StateOne const &rhs) const {
    return (species == rhs.species) && (n == rhs.n) && (l == rhs.l) && (j == rhs.j) && (m == rhs.m);
}
bool StateOne::operator^(StateOne const &rhs) const {
    return (species == rhs.species) && (rhs.n == ARB || n == rhs.n) &&
        (rhs.l == ARB || l == rhs.l) && (rhs.j == ARB || j == rhs.j) &&
        (rhs.m == ARB || m == rhs.m);
}
bool StateOne::operator!=(StateOne const &rhs) const {
    return (species != rhs.species) || (n != rhs.n) || (l != rhs.l) || (j != rhs.j) || (m != rhs.m);
}
bool StateOne::operator<(const StateOne &rhs) const {
    return (species < rhs.species) ||
        ((species == rhs.species) &&
         ((n < rhs.n) ||
          ((n == rhs.n) &&
           ((l < rhs.l) || ((l == rhs.l) && ((j < rhs.j) || ((j == rhs.j) && (m < rhs.m))))))));
}
bool StateOne::operator<=(const StateOne &rhs) const { return (*this < rhs) || (*this == rhs); }

// Utility methods
void StateOne::analyzeSpecies() {
    s = 0.5;
    element = species;
    if (std::isdigit(species.back()) != 0) {
        s = ((species.back() - '0') - 1) / 2.;
        element = species.substr(0, species.size() - 1);
    }
}
void StateOne::shouldBeArtificial(bool opinion) const {
    if (this->isArtificial() != opinion) {
        throw std::runtime_error("The state does not have this property.");
    }
}

////////////////////////////////////////////////////////////////////
/// Implementation of two-atom state ///////////////////////////////
////////////////////////////////////////////////////////////////////

StateTwo::StateTwo(std::array<std::string, 2> species, std::array<int, 2> n, std::array<int, 2> l,
                   std::array<float, 2> j, std::array<float, 2> m)
    : state_array({{StateOne(species[0], n[0], l[0], j[0], m[0]),
                    StateOne(species[1], n[1], l[1], j[1], m[1])}}) {
    hashvalue = 0;
    utils::hash_combine(hashvalue, state_array[0].getHash());
    utils::hash_combine(hashvalue, state_array[1].getHash());
}
StateTwo::StateTwo(std::array<std::string, 2> label)
    : state_array({{StateOne(label[0]), StateOne(label[1])}}) {
    hashvalue = 0;
    utils::hash_combine(hashvalue, state_array[0].getHash());
    utils::hash_combine(hashvalue, state_array[1].getHash());
}
StateTwo::StateTwo(StateOne first_state, StateOne second_state)
    : state_array({{std::move(first_state), std::move(second_state)}}) {
    hashvalue = 0;
    utils::hash_combine(hashvalue, state_array[0].getHash());
    utils::hash_combine(hashvalue, state_array[1].getHash());
}

// Methods for printing the state
std::ostream &operator<<(std::ostream &out, const StateTwo &state) {
    out << state.state_array[0] << state.state_array[1];
    return out;
}

std::string StateTwo::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

// Getters
std::array<int, 2> StateTwo::getN() const {
    return {{state_array[0].getN(), state_array[1].getN()}};
}
std::array<int, 2> StateTwo::getL() const {
    return {{state_array[0].getL(), state_array[1].getL()}};
}
std::array<float, 2> StateTwo::getJ() const {
    return {{state_array[0].getJ(), state_array[1].getJ()}};
}
std::array<float, 2> StateTwo::getM() const {
    return {{state_array[0].getM(), state_array[1].getM()}};
}
std::array<float, 2> StateTwo::getS() const {
    return {{state_array[0].getS(), state_array[1].getS()}};
}
std::array<std::string, 2> StateTwo::getSpecies() const {
    return {{state_array[0].getSpecies(), state_array[1].getSpecies()}};
}
std::array<std::string, 2> StateTwo::getElement() const {
    return {{state_array[0].getElement(), state_array[1].getElement()}};
}
double StateTwo::getEnergy() const {
    return state_array[0].getEnergy() + state_array[1].getEnergy();
}
double StateTwo::getEnergy(MatrixElementCache &cache) const {
    return state_array[0].getEnergy(cache) + state_array[1].getEnergy(cache);
}
std::array<double, 2> StateTwo::getNStar() const {
    return {{state_array[0].getNStar(), state_array[1].getNStar()}};
}
std::array<double, 2> StateTwo::getNStar(MatrixElementCache &cache) const {
    return {{state_array[0].getNStar(cache), state_array[1].getNStar(cache)}};
}
double StateTwo::getLeRoyRadius(MatrixElementCache &cache) const {
    return 2 *
        (std::sqrt(cache.getRadial(state_array[0], state_array[0], 2)) +
         std::sqrt(cache.getRadial(state_array[1], state_array[1], 2)));
}
std::array<std::string, 2> StateTwo::getLabel() const {
    return {{state_array[0].getLabel(), state_array[1].getLabel()}};
}
std::array<bool, 2> StateTwo::isArtificial() const {
    return {{state_array[0].isArtificial(), state_array[1].isArtificial()}};
}
std::array<bool, 2> StateTwo::isGeneralized() const {
    return {{state_array[0].isGeneralized(), state_array[1].isGeneralized()}};
}

const int &StateTwo::getN(int idx) const { return state_array[idx].getN(); }
const int &StateTwo::getL(int idx) const { return state_array[idx].getL(); }
const float &StateTwo::getJ(int idx) const { return state_array[idx].getJ(); }
const float &StateTwo::getM(int idx) const { return state_array[idx].getM(); }
const float &StateTwo::getS(int idx) const { return state_array[idx].getS(); }
const std::string &StateTwo::getSpecies(int idx) const { return state_array[idx].getSpecies(); }
const std::string &StateTwo::getElement(int idx) const { return state_array[idx].getElement(); }
double StateTwo::getEnergy(int idx) const { return state_array[idx].getEnergy(); }
double StateTwo::getEnergy(int idx, MatrixElementCache &cache) const {
    return state_array[idx].getEnergy(cache);
}
double StateTwo::getNStar(int idx) const { return state_array[idx].getNStar(); }
double StateTwo::getNStar(int idx, MatrixElementCache &cache) const {
    return state_array[idx].getNStar(cache);
}
const std::string &StateTwo::getLabel(int idx) const { return state_array[idx].getLabel(); }
bool StateTwo::isArtificial(int idx) const { return state_array[idx].isArtificial(); }
bool StateTwo::isGeneralized(int idx) const { return state_array[idx].isGeneralized(); }

const StateOne &StateTwo::getFirstState() const { return state_array[0]; }
const StateOne &StateTwo::getSecondState() const { return state_array[1]; }

const size_t &StateTwo::getHash() const { return hashvalue; }

StateTwo StateTwo::getReflected() const {
    return StateTwo(this->getSpecies(), this->getN(), this->getL(), this->getJ(),
                    {{-this->getM(0), -this->getM(1)}});
}

// Comparators
bool StateTwo::operator==(StateTwo const &rhs) const {
    return (state_array[0] == rhs.state_array[0]) && (state_array[1] == rhs.state_array[1]);
}
bool StateTwo::operator^(StateTwo const &rhs) const {
    return (state_array[0] ^ rhs.state_array[0]) && (state_array[1] ^ rhs.state_array[1]);
}
bool StateTwo::operator!=(StateTwo const &rhs) const {
    return (state_array[0] != rhs.state_array[0]) || (state_array[1] != rhs.state_array[1]);
}
bool StateTwo::operator<(const StateTwo &rhs) const {
    return (state_array[0] < rhs.state_array[0]) ||
        ((state_array[0] == rhs.state_array[0]) && (state_array[1] < rhs.state_array[1]));
}
bool StateTwo::operator<=(const StateTwo &rhs) const { return (*this < rhs) || (*this == rhs); }
