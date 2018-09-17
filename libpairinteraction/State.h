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
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/string.hpp>
#include <boost/variant.hpp>

////////////////////////////////////////////////////////////////////
/// Internally used state classes //////////////////////////////////
////////////////////////////////////////////////////////////////////

class StateInternalBase {
public:
    StateInternalBase() = default;
    friend std::ostream &operator<<(std::ostream &out, const StateInternalBase &state);
    virtual ~StateInternalBase() = default;
    virtual bool isArtificial() const = 0;
    bool operator==(StateInternalBase const &rhs) const;
    bool operator^(StateInternalBase const &rhs) const; // subset
    bool operator!=(StateInternalBase const &rhs) const;
    bool operator<(StateInternalBase const &rhs) const;

protected:
    virtual void printState(std::ostream &out) const = 0;
    virtual bool operatorEqual(StateInternalBase const &rhsb) const = 0;
    virtual bool operatorUnequal(StateInternalBase const &rhsb) const = 0;
    virtual bool operatorSubset(StateInternalBase const &rhsb) const = 0;
    virtual bool operatorLess(StateInternalBase const &rhsb) const = 0;

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive & /*ar*/, const unsigned int /*version*/) {}
};

class StateInternalFinestructure : public StateInternalBase {
public:
    StateInternalFinestructure() = default;
    StateInternalFinestructure(std::string species, int n, int l, float j, float m);
    bool isArtificial() const override;
    const std::string &getSpecies() const;
    const std::string &getElement() const;
    const int &getN() const;
    const int &getL() const;
    const float &getJ() const;
    const float &getM() const;
    const float &getS() const;
    double getEnergy() const;
    double getNStar() const;

protected:
    void printState(std::ostream &out) const override;
    bool operatorEqual(StateInternalBase const &rhsb) const override;
    bool operatorUnequal(StateInternalBase const &rhsb) const override;
    bool operatorSubset(StateInternalBase const &rhsb) const override;
    bool operatorLess(StateInternalBase const &rhsb) const override;

private:
    std::string species, element;
    int n, l;
    float j, m, s;

    void analyzeSpecies();

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/) {
        ar &boost::serialization::base_object<StateInternalBase>(*this);
        ar &species &n &l &j &m;
    }
};

class StateInternalArtificial : public StateInternalBase {
public:
    StateInternalArtificial() = default;
    StateInternalArtificial(std::string label);
    bool isArtificial() const override;
    const std::string &getLabel() const;

protected:
    void printState(std::ostream &out) const override;
    bool operatorEqual(StateInternalBase const &rhsb) const override;
    bool operatorUnequal(StateInternalBase const &rhsb) const override;
    bool operatorSubset(StateInternalBase const &rhsb) const override;
    bool operatorLess(StateInternalBase const &rhsb) const override;

private:
    std::string label;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/) {
        ar &boost::serialization::base_object<StateInternalBase>(*this);
        ar &label;
    }
};

// TODO make StateInternalBase etc. a nested class of StateOne

////////////////////////////////////////////////////////////////////
/// \brief One-atom state
/// This class implements a one-atom state. It can either be a
/// Rydberg state in the fine structure basis or an artificial state
/// specified by a label.
////////////////////////////////////////////////////////////////////

class StateOne {
public:
    StateOne() = default;
    StateOne(std::string species, int n, int l, float j, float m);
    StateOne(std::string label);

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

    // Comparators
    bool operator==(StateOne const &rhs) const;
    bool operator^(StateOne const &rhs) const; // subset
    bool operator!=(StateOne const &rhs) const;
    bool operator<(StateOne const &rhs) const;

private:
    std::shared_ptr<StateInternalBase>
        state_ptr; // TODO use unique_ptr (requires writing of copy constructor)

    // Method for serialization
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/) {
        ar &state_ptr;
    }

    // Utility methods
    template <class D>
    std::shared_ptr<D> castState(const std::shared_ptr<StateInternalBase> &p) const {
        auto d = std::dynamic_pointer_cast<D>(p);
        if (d == nullptr) {
            throw std::runtime_error("The state does not have this property.");
        }
        return d;
    }
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
    StateTwo(std::array<std::string, 2> species, std::array<int, 2> n, std::array<int, 2> l,
             std::array<float, 2> j, std::array<float, 2> m);
    StateTwo(std::array<std::string, 2> label);
    StateTwo(StateOne first_state, StateOne second_state);

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

    // Comparators
    bool operator==(StateTwo const &rhs) const;
    bool operator^(StateTwo const &rhs) const; // subset
    bool operator!=(StateTwo const &rhs) const;
    bool operator<(StateTwo const &rhs) const;

private:
    std::array<StateOne, 2> state_array;

    // Method for serialization
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int /*version*/) {
        ar &state_array;
    }
};

////////////////////////////////////////////////////////////////////
/// Hashers ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

#ifndef SWIG

namespace std {

template <>
struct hash<StateOne> {
    size_t operator()(const StateOne &s) const {
        if (s.isArtificial()) {
            return boost::hash_value(s.getLabel());
        }
        std::size_t seed = 0;
        boost::hash_combine(seed, s.getSpecies());
        boost::hash_combine(seed, s.getN());
        boost::hash_combine(seed, s.getL());
        boost::hash_combine(seed, s.getJ());
        boost::hash_combine(seed, s.getM());
        return seed;
    }
};

template <>
struct hash<StateTwo> {
    size_t operator()(const StateTwo &s) const {
        std::size_t seed = 0;
        for (int idx = 0; idx < 2; ++idx) {
            if (s.isArtificial(idx)) {
                boost::hash_combine(seed, s.getLabel(idx));
            } else {
                boost::hash_combine(seed, s.getSpecies(idx));
                boost::hash_combine(seed, s.getN(idx));
                boost::hash_combine(seed, s.getL(idx));
                boost::hash_combine(seed, s.getJ(idx));
                boost::hash_combine(seed, s.getM(idx));
            }
        }
        return seed;
    }
};

} // namespace std
#endif

#endif
