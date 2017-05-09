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

#ifndef STATE_H
#define STATE_H

#include "dtypes.h"
#include "utils.h"
#include "ConfParser.h"
#include "QuantumDefect.h"

#include <array>
#include <string>
#include <iostream>
#include <iomanip>

#include <boost/functional/hash.hpp>

class State {
public:
    State(idx_t idx) : idx(idx) { }
    idx_t idx;
};

class StateOne : public State {
public:
    StateOne(std::string element, int n, int l, float j, float m) : State(0), element(element), n(n), l(l), j(j), m(m) {
    }
    StateOne() : State(0), n(0), l(0), j(0), m(0) { }

    StateOne(idx_t idx, int n, int l, float s, float j, float m) : State(idx), n(n), l(l), s(s), j(j), m(m) {

    }
    StateOne(idx_t idx, int n, int l, float j, float m) : State(idx), n(n), l(l), s(0.5), j(j), m(m) {

    }
    StateOne(int n, int l, float j, float m) : State(0), n(n), l(l), s(0.5), j(j), m(m) {

    }

    StateOne(int n, int l, float s, float j, float m) : State(0), n(n), l(l), s(s), j(j), m(m) { }
    friend std::ostream& operator<< (std::ostream &out, const StateOne &state) {
        out << "n  =" << std::setw(3) << state.n << ",   ";
        out << "l  =" << std::setw(3) << state.l << ",   ";
        out << "j  =" << std::setprecision(2) << std::setw(5) << state.j << ",   ";
        out << "m  =" << std::setprecision(2) << std::setw(5) << state.m;
        return out;
    }

    friend bool operator== (const StateOne& s1, const StateOne& s2)
    {
        return ((s1.n == s2.n) && (s1.l == s2.l)  && (s1.s == s2.s)  && (s1.j == s2.j)  && (s1.m == s2.m)); // TODO compare element in addition
    }

    friend bool operator!= (const StateOne& s1, const StateOne& s2)
    {
        return ((s1.n != s2.n) || (s1.l != s2.l)  || (s1.s != s2.s)  || (s1.j != s2.j)  || (s1.m != s2.m));
    }

    friend bool operator< (const StateOne& s1, const StateOne& s2)
    {
        return (s1.idx  < s2.idx);
    }

    friend bool operator> (const StateOne& s1, const StateOne& s2)
    {
        return (s1.idx > s2.idx);
    }

    double getEnergy() const {
        return energy_level(element, n, l, j);
    }
    std::string getElement() const {
        return element;
    }
    int getN() const {
        return n;
    }
    int getL() const {
        return l;
    }
    float getJ() const {
        return j;
    }
    float getM() const {
        return m;
    }
    void setElement(std::string input){
        element = input;
    }
    void setN(int input){
        n = input;
    }
    void setL(int input){
        l = input;
    }
    void setJ(float input){
        j = input;
    }
    void setM(float input){
        m = input;
    }

    std::string element;
    int n, l;
    float s, j, m;
};

class StateTwo : public State {
public:
    StateTwo(std::array<std::string, 2> element, std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> j, std::array<float, 2> m) : State(0), element(element), n(n), l(l), j(j), m(m) {
    }
    StateTwo(const StateOne &s1, const StateOne &s2) : State(0), element({{s1.element, s2.element}}), n({{s1.n, s2.n}}), l({{s1.l, s2.l}}), j({{s1.j, s2.j}}), m({{s1.m, s2.m}}) {
    }
    StateOne getFirstState() const {
        return StateOne(element[0], n[0], l[0], j[0], m[0]);
    }
    StateOne getSecondState() const {
        return StateOne(element[1], n[1], l[1], j[1], m[1]);
    }
    void setFirstState(const StateOne &s) {
        element[0] = s.element;
        n[0] = s.n;
        l[0] = s.l;
        j[0] = s.j;
        m[0] = s.m;
    }
    void setSecondState(const StateOne &s) {
        element[1] = s.element;
        n[1] = s.n;
        l[1] = s.l;
        j[1] = s.j;
        m[1] = s.m;
    }

    StateTwo() : State(0), n({{0,0}}), l({{0,0}}), s({{0.5,0.5}}), j({{0,0}}), m({{0,0}}) { }
    StateTwo(idx_t idx, std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> s, std::array<float, 2> j, std::array<float, 2> m) : State(idx), n(n), l(l), s(s), j(j), m(m) { }
    StateTwo(std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> j, std::array<float, 2> m) : State(0), n(n), l(l), j(j), m(m) { }
    StateTwo(idx_t idx, const StateOne &a, const StateOne &b) : State(idx), n({{a.n,b.n}}), l({{a.l,b.l}}), s({{a.s,b.s}}), j({{a.j,b.j}}), m({{a.m,b.m}}) { }
    StateOne first() const {
        return StateOne(element[0], n[0], l[0], j[0], m[0]);
    }
    StateOne second() const {
        return StateOne(element[1], n[1], l[1], j[1], m[1]);
    }
    friend std::ostream& operator<< (std::ostream &out, const StateTwo &state) {
        for (size_t i = 0; i < 2; ++i) {
            out << "n" << i << " =" << std::setw(3) << state.n[i] << ",   ";
            out << "l" << i << " =" << std::setw(3) << state.l[i] << ",   ";
            out << "j" << i << " =" << std::setprecision(2) << std::setw(5) << state.j[i] << ",   ";
            out << "m" << i << " =" << std::setprecision(2) << std::setw(5) << state.m[i];
            if (i == 0) out << ",   ";
        }
        return out;
    }
    friend bool operator== (const StateTwo& s1, const StateTwo& s2)
    {
        return ((s1.n[0] == s2.n[0]) && (s1.l[0] == s2.l[0])  && (s1.s[0] == s2.s[0])  && (s1.j[0] == s2.j[0])  && (s1.m[0] == s2.m[0]) &&
                (s1.n[1] == s2.n[1]) && (s1.l[1] == s2.l[1])  && (s1.s[1] == s2.s[1])  && (s1.j[1] == s2.j[1])  && (s1.m[1] == s2.m[1]));
    }
    friend bool operator!= (const StateTwo& s1, const StateTwo& s2)
    {
        return ((s1.n[0] != s2.n[0]) || (s1.l[0] != s2.l[0])  || (s1.s[0] != s2.s[0])  || (s1.j[0] != s2.j[0])  || (s1.m[0] != s2.m[0]) ||
                (s1.n[1] != s2.n[1]) || (s1.l[1] != s2.l[1])  || (s1.s[1] != s2.s[1])  || (s1.j[1] != s2.j[1])  || (s1.m[1] != s2.m[1]));
    }
    /*StateOne first() const {
        return StateOne(idx, n[0], l[0], s[0], j[0], m[0]);
    }
    StateOne second() const {
        return StateOne(idx, n[1], l[1], s[1], j[1], m[1]);
    }*/
    StateTwo order() {
        if ((n[0] < n[1]) || ((n[0] == n[1]) &&
                              ((l[0] < l[1]) || ((l[0] == l[1]) &&
                                                 ((s[0] < s[1]) || ((s[0] == s[1]) &&
                                                                    ((j[0] < j[1]) || ((j[0] == j[1]) &&
                                                                                       (m[0] <= m[1]))))))))) {
            return *this;
        } else {
            return StateTwo(this->second(),this->first());
        }
    }
    double getEnergy() const {
        return this->first().getEnergy()+this->second().getEnergy();
    }
    std::array<std::string, 2> getElement() const {
        return element;
    }
    std::array<int, 2> getN() const {
        return n;
    }
    std::array<int, 2> getL() const {
        return l;
    }
    std::array<float, 2> getJ() const {
        return j;
    }
    std::array<float, 2> getM() const {
        return m;
    }
    void setElement(std::array<std::string, 2> input) {
        element = input;
    }
    void setN(std::array<int, 2> input) {
        n = input;
    }
    void setL(std::array<int, 2> input) {
        l = input;
    }
    void setJ(std::array<float, 2> input) {
        j = input;
    }
    void setM(std::array<float, 2> input) {
        m = input;
    }

    std::array<std::string, 2> element;
    std::array<int, 2> n, l;
    std::array<float, 2> s, j, m;
};

namespace std {
template <> struct hash<StateOne>
{
    size_t operator()(const StateOne & s) const
    {
        size_t seed = 0;
        boost::hash_combine(seed, s.n);
        boost::hash_combine(seed, s.l);
        boost::hash_combine(seed, s.s);
        boost::hash_combine(seed, s.j);
        boost::hash_combine(seed, s.m);
        return seed;
    }
};
template <> struct hash<StateTwo>
{
    size_t operator()(const StateTwo & s) const
    {
        size_t seed = 0;
        boost::hash_combine(seed, s.n[0]);
        boost::hash_combine(seed, s.l[0]);
        boost::hash_combine(seed, s.s[0]);
        boost::hash_combine(seed, s.j[0]);
        boost::hash_combine(seed, s.m[0]);
        boost::hash_combine(seed, s.n[1]);
        boost::hash_combine(seed, s.l[1]);
        boost::hash_combine(seed, s.s[1]);
        boost::hash_combine(seed, s.j[1]);
        boost::hash_combine(seed, s.m[1]);
        return seed;
    }
};
}

#endif
