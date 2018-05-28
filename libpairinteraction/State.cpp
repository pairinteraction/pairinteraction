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
#include "State.h"
#include "QuantumDefect.h"

#include <cctype>
#include <array>
#include <string>
#include <iostream>
#include <iomanip>

///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// Implementation of StateOne +++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

StateOne::StateOne(std::string element, int n, int l, float j, float m)
    : State(0), species(element), n(n), l(l), j(j), m(m)
{
    this->analyzeSpecies();
}

StateOne::StateOne()
    : State(0), species(""), n(0), l(0), j(0), m(0)
{
    this->analyzeSpecies();
}

StateOne::StateOne(idx_t idx, int n, int l, float j, float m)
    : State(idx), n(n), l(l), j(j), m(m)
{
    this->analyzeSpecies();
}

StateOne::StateOne(int n, int l, float j, float m)
    : State(0), n(n), l(l), j(j), m(m)
{
    this->analyzeSpecies();
}

std::ostream& operator<< (std::ostream &out, const StateOne &state)
{
    out << "n  =" << std::setw(3) << state.n << ",   ";
    out << "l  =" << std::setw(3) << state.l << ",   ";
    out << "j  =" << std::setprecision(4) << std::setw(5) << state.j << ",   ";
    out << "m  =" << std::setprecision(4) << std::setw(5) << state.m;
    return out;
}

bool StateOne::operator==(StateOne const& rhs) const
{
    // TODO use elements, too?
    return (n == rhs.n) && (l == rhs.l) && (j == rhs.j)  && (m == rhs.m);
}

bool StateOne::operator^(StateOne const& rhs) const{ // subset // TODO is there a better operator to use?
    return (rhs.n == ARB || n == rhs.n) && (rhs.l == ARB || l == rhs.l) && (rhs.j == ARB || j == rhs.j)  && (rhs.m == ARB || m == rhs.m);
}

bool StateOne::operator!=(StateOne const& rhs) const
{
    return ((n != rhs.n) || (l != rhs.l)  || (j != rhs.j)  || (m != rhs.m));
}

bool StateOne::operator<(const StateOne& rhs) const
{
    // TODO use elements, too?
    return ((n < rhs.n) || ((n == rhs.n) &&
                            ((l < rhs.l) || ((l == rhs.l) &&
                                             ((j < rhs.j) || ((j == rhs.j) &&
                                                              (m < rhs.m)))))));
}

bool StateOne::operator>(StateOne const& rhs) const // TODO remove this operator
{
    return (idx > rhs.idx);
}

double StateOne::getEnergy() const
{
    return energy_level(species, n, l, j);
}

double StateOne::getNStar() const
{
    return nstar(species, n, l, j);
}


////////////////////////////////////////////////////////////////////
/// Utility methods ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

void StateOne::analyzeSpecies() {
    s = 0.5;
    element = species;

    if (std::isdigit(species.back())) {
        s = (std::atoi(&species.back())-1)/2.;
        element = species.substr(0, species.size()-1);
    }
}

///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// Implementation of StateTwo +++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

StateTwo::StateTwo()
    : State(0), species({{"",""}}), n({{0,0}}), l({{0,0}}), j({{0,0}}), m({{0,0}})
{
    this->analyzeSpecies();
}

StateTwo::StateTwo(std::array<std::string, 2> element, std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> j, std::array<float, 2> m)
    : State(0), species(element), n(n), l(l), j(j), m(m)
{
    this->analyzeSpecies();
}

StateTwo::StateTwo(const StateOne &s1, const StateOne &s2)
    : State(0), species({{s1.species, s2.species}}), n({{s1.n, s2.n}}), l({{s1.l, s2.l}}), j({{s1.j, s2.j}}), m({{s1.m, s2.m}})
{
    this->analyzeSpecies();
}

StateTwo::StateTwo(idx_t idx, std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> j, std::array<float, 2> m)
    : State(idx), n(n), l(l), j(j), m(m)
{
    this->analyzeSpecies();
}

StateTwo::StateTwo(std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> j, std::array<float, 2> m)
    : State(0), n(n), l(l), j(j), m(m)
{
    this->analyzeSpecies();
}

StateTwo::StateTwo(idx_t idx, const StateOne &a, const StateOne &b)
    : State(idx), n({{a.n,b.n}}), l({{a.l,b.l}}), j({{a.j,b.j}}), m({{a.m,b.m}})
{
    this->analyzeSpecies();
}


StateOne StateTwo::getFirstState() const
{
    return StateOne(species[0], n[0], l[0], j[0], m[0]);
}
void StateTwo::setFirstState(StateOne const& s)
{
    species[0] = s.species;
    n[0] = s.n;
    l[0] = s.l;
    j[0] = s.j;
    m[0] = s.m;
}

StateOne StateTwo::getSecondState() const
{
    return StateOne(species[1], n[1], l[1], j[1], m[1]);
}
void StateTwo::setSecondState(StateOne const& s)
{
    species[1] = s.species;
    n[1] = s.n;
    l[1] = s.l;
    j[1] = s.j;
    m[1] = s.m;
}

StateOne StateTwo::first() const
{
    return StateOne(species[0], n[0], l[0], j[0], m[0]);
}
StateOne StateTwo::second() const
{
    return StateOne(species[1], n[1], l[1], j[1], m[1]);
}

std::ostream& operator<< (std::ostream &out, const StateTwo &state) {
    for (size_t i = 0; i < 2; ++i) {
        out << "n" << i << " =" << std::setw(3) << state.n[i] << ",   ";
        out << "l" << i << " =" << std::setw(3) << state.l[i] << ",   ";
        out << "j" << i << " =" << std::setprecision(4) << std::setw(5) << state.j[i] << ",   ";
        out << "m" << i << " =" << std::setprecision(4) << std::setw(5) << state.m[i];
        if (i == 0) out << ",   ";
    }
    return out;
}

bool StateTwo::operator==(const StateTwo& rhs) const
{
    return (n[0] == rhs.n[0]) && (l[0] == rhs.l[0])  && (j[0] == rhs.j[0])  && (m[0] == rhs.m[0]) &&
            (n[1] == rhs.n[1]) && (l[1] == rhs.l[1])  && (j[1] == rhs.j[1])  && (m[1] == rhs.m[1]);
}

bool StateTwo::operator^(const StateTwo& rhs) const{ // subset // TODO is there a better operator to use?
    return (rhs.n[0] == ARB || n[0] == rhs.n[0]) && (rhs.l[0] == ARB || l[0] == rhs.l[0])  && (rhs.j[0] == ARB || j[0] == rhs.j[0])  && (rhs.m[0] == ARB || m[0] == rhs.m[0]) &&
            (rhs.n[1] == ARB || n[1] == rhs.n[1]) && (rhs.l[1] == ARB || l[1] == rhs.l[1])  && (rhs.j[1] == ARB || j[1] == rhs.j[1])  && (rhs.m[1] == ARB || m[1] == rhs.m[1]);
}

bool StateTwo::operator!=(const StateTwo& rhs) const
{
    return (n[0] != rhs.n[0]) || (l[0] != rhs.l[0]) || (j[0] != rhs.j[0])  || (m[0] != rhs.m[0]) ||
            (n[1] != rhs.n[1]) || (l[1] != rhs.l[1]) || (j[1] != rhs.j[1])  || (m[1] != rhs.m[1]);
}

bool StateTwo::operator<(const StateTwo& rhs) const
{
    return ((this->first() < rhs.first()) || ((this->first() == rhs.first()) &&
                                              (this->second() < rhs.second())));
}

StateTwo StateTwo::order() { // TODO use element, too?
    if ((n[0] < n[1]) || ((n[0] == n[1]) &&
                          ((l[0] < l[1]) || ((l[0] == l[1]) &&
                                             ((j[0] < j[1]) || ((j[0] == j[1]) &&
                                                                (m[0] <= m[1]))))))) {
        return *this;
    } else {
        return StateTwo(this->second(),this->first());
    }
}

double StateTwo::getEnergy() const {
    return this->first().getEnergy()+this->second().getEnergy();
}

std::array<double, 2> StateTwo::getNStar() const {
    return {{this->first().getNStar(), this->second().getNStar()}};
}

////////////////////////////////////////////////////////////////////
/// Utility methods ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

void StateTwo::analyzeSpecies() {
    for (size_t i = 0; i < 2; ++i) {
        s[i] = 0.5;
        element[i] = species[i];

        if (std::isdigit(species[i].back())) {
            s[i] = (std::atoi(&species[i].back())-1)/2.;
            element[i] = species[i].substr(0, species[i].size()-1);
        }
    }
}
