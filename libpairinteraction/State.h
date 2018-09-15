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

#include <boost/serialization/array.hpp>
#include <boost/serialization/string.hpp>
#include <boost/variant.hpp>

/** \brief %Base class for states
 *
 * This class is the base class for states specified in the fine structure basis.
 */
class State {
public:
    State(idx_t idx) : idx(idx) {}
    idx_t idx;
};

/** \brief %One-atom Rydberg state
 *
 * This class implements a one-atom Rydberg state.
 */
class StateOne : public State {
public:
    // These are public to allow direct access.  This violates the
    // open/closed principle and is a sign of code smell.
    std::string species, element;
    int n{0};
    int l{0};
    float j{0};
    float m{0};
    float s;

    StateOne();
    StateOne(std::string element, int n, int l, float j, float m);
    StateOne(idx_t idx, int n, int l, float j, float m);
    StateOne(int n, int l, float j, float m);

    friend std::ostream &operator<<(std::ostream &out, const StateOne &state);

    bool operator==(StateOne const & /*rhs*/) const;
    bool operator^(StateOne const & /*rhs*/) const; // subset
    bool operator!=(StateOne const & /*rhs*/) const;
    bool operator<(StateOne const & /*rhs*/) const;

    double getEnergy() const;
    double getNStar() const;

    std::string getSpecies() const;
    int getN() const;
    int getL() const;
    float getJ() const;
    float getM() const;

private:
    ////////////////////////////////////////////////////////////////////
    /// Utility methods ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    void analyzeSpecies();

    ////////////////////////////////////////////////////////////////////
    /// Method for serialization ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        (void)version;

        ar &species &element &s &n &l &j &m;
    }
};

/** \brief %Artificial one-atom Rydberg state
 *
 * This class implements an artificial one-atom Rydberg state.
 */
class StateOneArtificial : private StateOne { // TODO private StateOne, StateArtificial
public:
    StateOneArtificial();
    StateOneArtificial(std::string label);
    friend std::ostream &operator<<(std::ostream &out, const StateOneArtificial &state);
    bool operator==(StateOneArtificial const &rhs) const;
    bool operator!=(StateOneArtificial const &rhs) const;
    bool operator<(const StateOneArtificial &rhs) const;
    std::string getLabel() const;

private:
    std::string label;

    // Method for serialization
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/) {
        ar &label;
    }
};

/** \brief %Two-atom Rydberg state
 *
 * This class implements a two-atom Rydberg state.
 */
class StateTwo
    : public State { // TODO define getters and setters, save a pair state as two single atom states
public:
    // These are public to allow direct access.  This violates the
    // open/closed principle and is a sign of code smell.
    std::array<std::string, 2> species, element;
    std::array<int, 2> n, l;
    std::array<float, 2> j, m, s;

    StateTwo();
    StateTwo(std::array<std::string, 2> element, std::array<int, 2> n, std::array<int, 2> l,
             std::array<float, 2> j, std::array<float, 2> m);
    StateTwo(const StateOne &s1, const StateOne &s2);
    StateTwo(idx_t idx, std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> j,
             std::array<float, 2> m);
    StateTwo(std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> j,
             std::array<float, 2> m);
    StateTwo(idx_t idx, const StateOne &a, const StateOne &b);

    StateOne getFirstState() const;
    StateOne getSecondState() const;
    void setFirstState(StateOne const & /*s*/);
    void setSecondState(StateOne const & /*s*/);

    StateOne first() const;
    StateOne second() const;

    friend std::ostream &operator<<(std::ostream &out, const StateTwo &state);

    bool operator==(StateTwo const & /*rhs*/) const;
    bool operator^(StateTwo const & /*rhs*/) const; // subset
    bool operator!=(StateTwo const & /*rhs*/) const;
    bool operator<(StateTwo const & /*rhs*/) const;

    StateTwo order();

    double getEnergy() const;
    std::array<double, 2> getNStar() const;

    std::array<std::string, 2> getSpecies() const;
    std::array<int, 2> getN() const;
    std::array<int, 2> getL() const;
    std::array<float, 2> getJ() const;
    std::array<float, 2> getM() const;

private:
    ////////////////////////////////////////////////////////////////////
    /// Utility methods ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    void analyzeSpecies();

    ////////////////////////////////////////////////////////////////////
    /// Method for serialization ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        (void)version;

        ar &species &element &s &n &l &j &m;
    }
};

/** \brief %Artificial two-atom Rydberg state
 *
 * This class implements an artificial two-atom Rydberg state.
 */
class StateTwoArtificial : private StateTwo {
public:
    StateTwoArtificial();
    StateTwoArtificial(std::array<std::string, 2> label);
    friend std::ostream &operator<<(std::ostream &out, const StateTwoArtificial &state);
    bool operator==(StateTwoArtificial const &rhs) const;
    bool operator!=(StateTwoArtificial const &rhs) const;
    bool operator<(const StateTwoArtificial &rhs) const;
    std::array<std::string, 2> getLabel() const;
    StateOneArtificial getFirstState() const;
    StateOneArtificial getSecondState() const;

private:
    std::array<std::string, 2> label;

    // Method for serialization
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/) {
        ar &label;
    }
};

#ifndef SWIG
namespace std {

template <>
struct hash<StateOne> {
    size_t operator()(const StateOne &s) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, s.n);
        boost::hash_combine(seed, s.l);
        boost::hash_combine(seed, s.j);
        boost::hash_combine(seed, s.m);
        return seed;
    }
};

template <>
struct hash<StateOneArtificial> {
    size_t operator()(const StateOneArtificial &s) const { return boost::hash_value(s.getLabel()); }
};

template <>
struct hash<StateTwo> {
    size_t operator()(const StateTwo &s) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, s.n);
        boost::hash_combine(seed, s.l);
        boost::hash_combine(seed, s.j);
        boost::hash_combine(seed, s.m);
        return seed;
    }
};

template <>
struct hash<StateTwoArtificial> {
    size_t operator()(const StateTwoArtificial &s) const { return boost::hash_value(s.getLabel()); }
};

} // namespace std
#endif

#endif
