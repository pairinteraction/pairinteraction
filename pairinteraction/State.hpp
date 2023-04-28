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

#ifndef STATE_H
#define STATE_H

#include "MatrixElementCache.hpp"
#include "dtypes.hpp"
#include "utils.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

// clang-format off
#if __has_include (<boost/serialization/version.hpp>)
#    include <boost/serialization/version.hpp>
#endif
#if __has_include (<boost/serialization/library_version_type.hpp>)
#    include <boost/serialization/library_version_type.hpp>
#endif
// clang-format on
#include <boost/serialization/array.hpp>
#include <boost/serialization/string.hpp>

class MatrixElementCache;

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

    // Methods for printing the state
    friend std::ostream &operator<<(std::ostream &out, const StateOne &state);
    std::string str() const;

    // Getters
    const int &getN() const;
    const int &getL() const;
    const float &getJ() const;
    const float &getM() const;
    const float &getS() const;
    const std::string &getSpecies() const;
    const std::string &getElement() const;
    double getEnergy() const;
    double getEnergy(MatrixElementCache &cache) const;
    double getNStar() const;
    double getNStar(MatrixElementCache &cache) const;
    const std::string &getLabel() const;
    bool isArtificial() const;
    bool isGeneralized() const;

    const size_t &getHash() const;

    StateOne getReflected() const;

    // Comparators
    bool operator==(StateOne const &rhs) const;
    bool operator^(StateOne const &rhs) const; // subset
    bool operator!=(StateOne const &rhs) const;
    bool operator<(StateOne const &rhs) const;
    bool operator<=(StateOne const &rhs) const;

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

    // Methods for printing the state
    friend std::ostream &operator<<(std::ostream &out, const StateTwo &state);
    std::string str() const;

    // Getters
    std::array<int, 2> getN() const;
    std::array<int, 2> getL() const;
    std::array<float, 2> getJ() const;
    std::array<float, 2> getM() const;
    std::array<float, 2> getS() const;
    std::array<std::string, 2> getSpecies() const;
    std::array<std::string, 2> getElement() const;
    double getEnergy() const;
    double getEnergy(MatrixElementCache &cache) const;
    std::array<double, 2> getNStar() const;
    std::array<double, 2> getNStar(MatrixElementCache &cache) const;
    double getLeRoyRadius(MatrixElementCache &cache) const;
    std::array<std::string, 2> getLabel() const;
    std::array<bool, 2> isArtificial() const;
    std::array<bool, 2> isGeneralized() const;

    const int &getN(int idx) const;
    const int &getL(int idx) const;
    const float &getJ(int idx) const;
    const float &getM(int idx) const;
    const float &getS(int idx) const;
    const std::string &getSpecies(int idx) const;
    const std::string &getElement(int idx) const;
    double getEnergy(int idx) const;
    double getEnergy(int idx, MatrixElementCache &cache) const;
    double getNStar(int idx) const;
    double getNStar(int idx, MatrixElementCache &cache) const;
    const std::string &getLabel(int idx) const;
    bool isArtificial(int idx) const;
    bool isGeneralized(int idx) const;

    const StateOne &getFirstState() const;
    const StateOne &getSecondState() const;

    const size_t &getHash() const;

    StateTwo getReflected() const;

    // Comparators
    bool operator==(StateTwo const &rhs) const;
    bool operator^(StateTwo const &rhs) const; // subset
    bool operator!=(StateTwo const &rhs) const;
    bool operator<(StateTwo const &rhs) const;
    bool operator<=(StateTwo const &rhs) const;

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
