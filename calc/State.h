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
#include <string>
#include <iostream>
#include <boost/functional/hash.hpp>

#include <boost/serialization/array.hpp>
#include <boost/serialization/string.hpp>

/** \brief %Base class for states
 *
 * This class is the base class for states specified in the fine structure basis.
 */ // TODO [dummystate]
class State {
public:
    State(idx_t idx) : idx(idx) { }
    idx_t idx;
};


/** \brief %One-atom Rydberg state
 *
 * This class implements a one-atom Rydberg state.
 */ // TODO [dummystate]
class StateOne : public State {
public:
    // These are public to allow direct access.  This violates the
    // open/closed principle and is a sign of code smell.
    std::wstring element;
    int n, l;
    float j, m;


    StateOne();
    StateOne(std::wstring element, int n, int l, float j, float m);
    StateOne(idx_t idx, int n, int l, float j, float m);
    StateOne(int n, int l, float j, float m);

    friend std::ostream& operator<<(std::ostream &out, const StateOne &state);

    bool operator==(StateOne const&) const;
    bool operator!=(StateOne const&) const;
    bool operator< (StateOne const&) const;
    bool operator> (StateOne const&) const;

    double getEnergy() const;
private:
    ////////////////////////////////////////////////////////////////////
    /// Method for serialization ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        (void)version;

        ar & element & n & l & j & m;
    }
};


#ifndef SWIG
namespace std {

    template <>
    struct hash<StateOne>
    {
        size_t operator()(const StateOne & s) const
        {
            size_t seed = 0; // TODO use element, too
            boost::hash_combine(seed, s.n);
            boost::hash_combine(seed, s.l);
            boost::hash_combine(seed, s.j);
            boost::hash_combine(seed, s.m);
            return seed;
        }
    };

}
#endif

/** \brief %Two-atom Rydberg state
 *
 * This class implements a two-atom Rydberg state.
 */ // TODO [dummystate]
class StateTwo : public State {
public:
    // These are public to allow direct access.  This violates the
    // open/closed principle and is a sign of code smell.
    std::array<std::wstring, 2> element;
    std::array<int, 2> n, l;
    std::array<float, 2> j, m;

    StateTwo();
    StateTwo(std::array<std::wstring, 2> element, std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> j, std::array<float, 2> m);
    StateTwo(const StateOne &s1, const StateOne &s2);
    StateTwo(idx_t idx, std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> j, std::array<float, 2> m);
    StateTwo(std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> j, std::array<float, 2> m);
    StateTwo(idx_t idx, const StateOne &a, const StateOne &b);

    StateOne getFirstState() const;
    StateOne getSecondState() const;
    void setFirstState(StateOne const&);
    void setSecondState(StateOne const&);

    StateOne first() const;
    StateOne second() const;

    friend std::ostream& operator<<(std::ostream &out, const StateTwo &state);

    bool operator==(StateTwo const&) const;
    bool operator!=(StateTwo const&) const;
    bool operator< (StateTwo const&) const;

    StateTwo order();

    double getEnergy() const;
private:
    ////////////////////////////////////////////////////////////////////
    /// Method for serialization ///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        (void)version;

        ar & element & n & l & j & m;
    }
};


#ifndef SWIG
namespace std {

    template <>
    struct hash<StateTwo>
    {
        size_t operator()(const StateTwo & s) const
        {
            size_t seed = 0;  // TODO use element, too
            boost::hash_combine(seed, s.n[0]);  // TODO do hash s.n etc. at once --> template for StateOne and StateTwo will work
            boost::hash_combine(seed, s.l[0]);
            boost::hash_combine(seed, s.j[0]);
            boost::hash_combine(seed, s.m[0]);
            boost::hash_combine(seed, s.n[1]);
            boost::hash_combine(seed, s.l[1]);
            boost::hash_combine(seed, s.j[1]);
            boost::hash_combine(seed, s.m[1]);
            return seed;
        }
    };

}
#endif


#endif
