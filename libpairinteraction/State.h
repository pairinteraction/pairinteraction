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

#include <array>
#include <boost/functional/hash.hpp>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <boost/serialization/array.hpp>
#include <boost/serialization/string.hpp>
#include <boost/variant.hpp>

////////////////////////////////////////////////////////////////////
/// \brief One-atom state
/// This class implements a one-atom state. It can either be a
/// Rydberg state in the fine structure basis or an artificial state
/// specified by a label.
////////////////////////////////////////////////////////////////////

class StateOne {
public:
    StateOne() = default;
    explicit StateOne(std::string species, int n, int l, float j, float m);
    explicit StateOne(std::string label);

    // Method for printing the state
    friend std::ostream &operator<<(std::ostream &out, const StateOne &state);

    // Getters
    const int &getN() const;
    const int &getL() const;
    const float &getJ() const;
    const float &getM() const;
    const float &getS() const;
    const std::string &getSpecies() const;
    const std::string &getElement() const;
    double getEnergy() const;
    double getNStar() const;
    const std::string &getLabel() const;
    bool isArtificial() const;

    const size_t &getHash() const;

    // Comparators
    bool operator==(StateOne const &rhs) const;
    bool operator^(StateOne const &rhs) const; // subset
    bool operator!=(StateOne const &rhs) const;
    bool operator<(StateOne const &rhs) const;

private:
    // TODO make the variables constant (requires load_construct_data, see
    // https://stackoverflow.com/questions/50603180/serialization-of-class-with-const-members-using-boost)
    std::string species, element;
    int n, l;
    float j, m, s;
    size_t hashvalue;

    // Method for serialization
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/) {
        ar &species &element &n &l &j &m &s &hashvalue;
    }

    // Utility methods
    void analyzeSpecies();
    void shouldBeArtificial(bool opinion) const;
};

////////////////////////////////////////////////////////////////////
/// \brief Two-atom state
/// This class implements a two-atom state. It can either be a
/// Rydberg state in the fine structure basis or an artificial state
/// specified by a label.
////////////////////////////////////////////////////////////////////

class StateTwo {
public:
    StateTwo() = default;
    explicit StateTwo(std::array<std::string, 2> species, std::array<int, 2> n,
                      std::array<int, 2> l, std::array<float, 2> j, std::array<float, 2> m);
    explicit StateTwo(std::array<std::string, 2> label); // TODO use &&label?
    explicit StateTwo(StateOne first_state, StateOne second_state);

    // Method for printing the state
    friend std::ostream &operator<<(std::ostream &out, const StateTwo &state);

    // Getters
    std::array<int, 2> getN() const;
    std::array<int, 2> getL() const;
    std::array<float, 2> getJ() const;
    std::array<float, 2> getM() const;
    std::array<float, 2> getS() const;
    std::array<std::string, 2> getSpecies() const;
    std::array<std::string, 2> getElement() const;
    double getEnergy() const;
    std::array<double, 2> getNStar() const;
    std::array<std::string, 2> getLabel() const;
    std::array<bool, 2> isArtificial() const;

    const int &getN(int idx) const;
    const int &getL(int idx) const;
    const float &getJ(int idx) const;
    const float &getM(int idx) const;
    const float &getS(int idx) const;
    const std::string &getSpecies(int idx) const;
    const std::string &getElement(int idx) const;
    double getEnergy(int idx) const;
    double getNStar(int idx) const;
    const std::string &getLabel(int idx) const;
    bool isArtificial(int idx) const;

    const StateOne &getFirstState() const;
    const StateOne &getSecondState() const;

    const size_t &getHash() const;

    // Comparators
    bool operator==(StateTwo const &rhs) const;
    bool operator^(StateTwo const &rhs) const; // subset
    bool operator!=(StateTwo const &rhs) const;
    bool operator<(StateTwo const &rhs) const;

private:
    // TODO make the variables constant (requires load_construct_data, see
    // https://stackoverflow.com/questions/50603180/serialization-of-class-with-const-members-using-boost)
    std::array<StateOne, 2> state_array;
    size_t hashvalue;

    // Method for serialization
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/) {
        ar &state_array &hashvalue;
    }
};

////////////////////////////////////////////////////////////////////
/// Hashers ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

#ifndef SWIG

namespace std {

template <>
struct hash<StateOne> {
    size_t operator()(const StateOne &s) const { return s.getHash(); }
};

template <>
struct hash<StateTwo> {
    size_t operator()(const StateTwo &s) const { return s.getHash(); }
};

} // namespace std
#endif

#endif
