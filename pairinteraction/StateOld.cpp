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

#include "StateOld.hpp"
#include "QuantumDefect.hpp"
#include "dtypes.hpp"

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

///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// Implementation of StateOne +++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

StateOneOld::StateOneOld(std::string element, int n, int l, float j, float m)
    : StateOld(0), species(std::move(element)), n(n), l(l), j(j), m(m) {
    this->analyzeSpecies();
}

StateOneOld::StateOneOld() : StateOld(0), species("") { this->analyzeSpecies(); }

StateOneOld::StateOneOld(idx_t idx, int n, int l, float j, float m)
    : StateOld(idx), n(n), l(l), j(j), m(m) {
    this->analyzeSpecies();
}

StateOneOld::StateOneOld(int n, int l, float j, float m) : StateOld(0), n(n), l(l), j(j), m(m) {
    this->analyzeSpecies();
}

std::ostream &operator<<(std::ostream &out, const StateOneOld &state) {
    out << "|" << state.species << ", ";

    out << state.n << " ";

    out << getMomentumLabel(state.l) << "_";

    if (std::ceil(state.j) == state.j) {
        out << state.j << ", ";
        out << "mj=" << state.m << ">";
    } else {
        out << 2 * state.j << "/2, ";
        out << "mj=" << 2 * state.m << "/2>";
    }

    return out;
}

bool StateOneOld::operator==(StateOneOld const &rhs) const {
    // TODO use elements, too?
    return (n == rhs.n) && (l == rhs.l) && (j == rhs.j) && (m == rhs.m);
}

bool StateOneOld::operator^(
    StateOneOld const &rhs) const { // subset // TODO is there a better operator to use?
    return (rhs.n == ARB || n == rhs.n) && (rhs.l == ARB || l == rhs.l) &&
        (rhs.j == ARB || j == rhs.j) && (rhs.m == ARB || m == rhs.m);
}

bool StateOneOld::operator!=(StateOneOld const &rhs) const {
    return ((n != rhs.n) || (l != rhs.l) || (j != rhs.j) || (m != rhs.m));
}

bool StateOneOld::operator<(const StateOneOld &rhs) const {
    // TODO use elements, too?
    return ((n < rhs.n) ||
            ((n == rhs.n) &&
             ((l < rhs.l) || ((l == rhs.l) && ((j < rhs.j) || ((j == rhs.j) && (m < rhs.m)))))));
}

double StateOneOld::getEnergy() const { return energy_level(species, n, l, j); }

double StateOneOld::getNStar() const { return nstar(species, n, l, j); }

std::string StateOneOld::getSpecies() const { return species; }

int StateOneOld::getN() const { return n; }

int StateOneOld::getL() const { return l; }

float StateOneOld::getJ() const { return j; }

float StateOneOld::getM() const { return m; }

////////////////////////////////////////////////////////////////////
/// Utility methods ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

void StateOneOld::analyzeSpecies() {
    s = 0.5;
    element = species;

    if (!species.empty() && std::isdigit(species.back()) != 0) {
        s = ((species.back() - '0') - 1) / 2.;
        element = species.substr(0, species.size() - 1);
    }
}

///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// Implementation of StateTwo +++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

StateTwoOld::StateTwoOld()
    : StateOld(0), species({{"", ""}}), n({{0, 0}}), l({{0, 0}}), j({{0, 0}}), m({{0, 0}}) {
    this->analyzeSpecies();
}

StateTwoOld::StateTwoOld(std::array<std::string, 2> element, std::array<int, 2> n,
                         std::array<int, 2> l, std::array<float, 2> j, std::array<float, 2> m)
    : StateOld(0), species(std::move(element)), n(n), l(l), j(j), m(m) {
    this->analyzeSpecies();
}

StateTwoOld::StateTwoOld(const StateOneOld &s1, const StateOneOld &s2)
    : StateOld(0), species({{s1.species, s2.species}}), n({{s1.n, s2.n}}), l({{s1.l, s2.l}}),
      j({{s1.j, s2.j}}), m({{s1.m, s2.m}}) {
    this->analyzeSpecies();
}

StateTwoOld::StateTwoOld(idx_t idx, std::array<int, 2> n, std::array<int, 2> l,
                         std::array<float, 2> j, std::array<float, 2> m)
    : StateOld(idx), n(n), l(l), j(j), m(m) {
    this->analyzeSpecies();
}

StateTwoOld::StateTwoOld(std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> j,
                         std::array<float, 2> m)
    : StateOld(0), n(n), l(l), j(j), m(m) {
    this->analyzeSpecies();
}

StateTwoOld::StateTwoOld(idx_t idx, const StateOneOld &a, const StateOneOld &b)
    : StateOld(idx), n({{a.n, b.n}}), l({{a.l, b.l}}), j({{a.j, b.j}}), m({{a.m, b.m}}) {
    this->analyzeSpecies();
}

StateOneOld StateTwoOld::getFirstState() const {
    return StateOneOld(species[0], n[0], l[0], j[0], m[0]);
}
void StateTwoOld::setFirstState(StateOneOld const &s) {
    species[0] = s.species;
    n[0] = s.n;
    l[0] = s.l;
    j[0] = s.j;
    m[0] = s.m;
}

StateOneOld StateTwoOld::getSecondState() const {
    return StateOneOld(species[1], n[1], l[1], j[1], m[1]);
}
void StateTwoOld::setSecondState(StateOneOld const &s) {
    species[1] = s.species;
    n[1] = s.n;
    l[1] = s.l;
    j[1] = s.j;
    m[1] = s.m;
}

StateOneOld StateTwoOld::first() const { return StateOneOld(species[0], n[0], l[0], j[0], m[0]); }
StateOneOld StateTwoOld::second() const { return StateOneOld(species[1], n[1], l[1], j[1], m[1]); }

std::ostream &operator<<(std::ostream &out, const StateTwoOld &state) {
    out << "|";
    for (size_t i = 0; i < 2; ++i) {
        out << state.species[i] << ", ";

        out << state.n[i] << " ";

        out << getMomentumLabel(state.l[i]) << "_";

        if (std::ceil(state.j[i]) == state.j[i]) {
            out << state.j[i] << ", ";
            out << "mj=" << state.m[i] << ">";
        } else {
            out << 2 * state.j[i] << "/2, ";
            out << "mj=" << 2 * state.m[i] << "/2";
        }

        if (i == 0) {
            out << "; ";
        }
    }
    out << ">";
    return out;
}

bool StateTwoOld::operator==(const StateTwoOld &rhs) const {
    return (n[0] == rhs.n[0]) && (l[0] == rhs.l[0]) && (j[0] == rhs.j[0]) && (m[0] == rhs.m[0]) &&
        (n[1] == rhs.n[1]) && (l[1] == rhs.l[1]) && (j[1] == rhs.j[1]) && (m[1] == rhs.m[1]);
}

bool StateTwoOld::operator^(
    const StateTwoOld &rhs) const { // subset // TODO is there a better operator to use?
    return (rhs.n[0] == ARB || n[0] == rhs.n[0]) && (rhs.l[0] == ARB || l[0] == rhs.l[0]) &&
        (rhs.j[0] == ARB || j[0] == rhs.j[0]) && (rhs.m[0] == ARB || m[0] == rhs.m[0]) &&
        (rhs.n[1] == ARB || n[1] == rhs.n[1]) && (rhs.l[1] == ARB || l[1] == rhs.l[1]) &&
        (rhs.j[1] == ARB || j[1] == rhs.j[1]) && (rhs.m[1] == ARB || m[1] == rhs.m[1]);
}

bool StateTwoOld::operator!=(const StateTwoOld &rhs) const {
    return (n[0] != rhs.n[0]) || (l[0] != rhs.l[0]) || (j[0] != rhs.j[0]) || (m[0] != rhs.m[0]) ||
        (n[1] != rhs.n[1]) || (l[1] != rhs.l[1]) || (j[1] != rhs.j[1]) || (m[1] != rhs.m[1]);
}

bool StateTwoOld::operator<(const StateTwoOld &rhs) const {
    return ((this->first() < rhs.first()) ||
            ((this->first() == rhs.first()) && (this->second() < rhs.second())));
}

StateTwoOld StateTwoOld::order() { // TODO use element, too?
    if ((n[0] < n[1]) ||
        ((n[0] == n[1]) &&
         ((l[0] < l[1]) ||
          ((l[0] == l[1]) && ((j[0] < j[1]) || ((j[0] == j[1]) && (m[0] <= m[1]))))))) {
        return *this;
    }
    return StateTwoOld(this->second(), this->first());
}

double StateTwoOld::getEnergy() const {
    return this->first().getEnergy() + this->second().getEnergy();
}

std::array<double, 2> StateTwoOld::getNStar() const {
    return {{this->first().getNStar(), this->second().getNStar()}};
}

std::array<std::string, 2> StateTwoOld::getSpecies() const { return species; }

std::array<int, 2> StateTwoOld::getN() const { return n; }

std::array<int, 2> StateTwoOld::getL() const { return l; }

std::array<float, 2> StateTwoOld::getJ() const { return j; }

std::array<float, 2> StateTwoOld::getM() const { return m; }

////////////////////////////////////////////////////////////////////
/// Utility methods ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

void StateTwoOld::analyzeSpecies() {
    for (size_t i = 0; i < 2; ++i) {
        s[i] = 0.5;
        element[i] = species[i];

        if (!species[i].empty() && std::isdigit(species[i].back()) != 0) {
            s[i] = ((species[i].back() - '0') - 1) / 2.;
            element[i] = species[i].substr(0, species[i].size() - 1);
        }
    }
}
