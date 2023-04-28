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

#ifndef STATEOLD_H
#define STATEOLD_H

#include "dtypes.hpp"
#include "utils.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <string>

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

/** \brief %Base class for states
 *
 * This class is the base class for states specified in the fine structure basis.
 */
class StateOld {
public:
    StateOld(idx_t idx) : idx(idx) {}
    idx_t idx;
};

/** \brief %One-atom Rydberg state
 *
 * This class implements a one-atom Rydberg state.
 */
class StateOneOld : public StateOld {
public:
    // These are public to allow direct access.  This violates the
    // open/closed principle and is a sign of code smell.
    std::string species, element;
    int n{0};
    int l{0};
    float j{0};
    float m{0};
    float s;

    StateOneOld();
    StateOneOld(std::string element, int n, int l, float j, float m);
    StateOneOld(idx_t idx, int n, int l, float j, float m);
    StateOneOld(int n, int l, float j, float m);

    friend std::ostream &operator<<(std::ostream &out, const StateOneOld &state);

    bool operator==(StateOneOld const & /*rhs*/) const;
    bool operator^(StateOneOld const & /*rhs*/) const; // subset
    bool operator!=(StateOneOld const & /*rhs*/) const;
    bool operator<(StateOneOld const & /*rhs*/) const;

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

/** \brief %Two-atom Rydberg state
 *
 * This class implements a two-atom Rydberg state.
 */
class StateTwoOld : public StateOld { // TODO define getters and setters, save a pair state as two
                                      // single atom states
public:
    // These are public to allow direct access.  This violates the
    // open/closed principle and is a sign of code smell.
    std::array<std::string, 2> species, element;
    std::array<int, 2> n, l;
    std::array<float, 2> j, m, s;

    StateTwoOld();
    StateTwoOld(std::array<std::string, 2> element, std::array<int, 2> n, std::array<int, 2> l,
                std::array<float, 2> j, std::array<float, 2> m);
    StateTwoOld(const StateOneOld &s1, const StateOneOld &s2);
    StateTwoOld(idx_t idx, std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> j,
                std::array<float, 2> m);
    StateTwoOld(std::array<int, 2> n, std::array<int, 2> l, std::array<float, 2> j,
                std::array<float, 2> m);
    StateTwoOld(idx_t idx, const StateOneOld &a, const StateOneOld &b);

    StateOneOld getFirstState() const;
    StateOneOld getSecondState() const;
    void setFirstState(StateOneOld const & /*s*/);
    void setSecondState(StateOneOld const & /*s*/);

    StateOneOld first() const;
    StateOneOld second() const;

    friend std::ostream &operator<<(std::ostream &out, const StateTwoOld &state);

    bool operator==(StateTwoOld const & /*rhs*/) const;
    bool operator^(StateTwoOld const & /*rhs*/) const; // subset
    bool operator!=(StateTwoOld const & /*rhs*/) const;
    bool operator<(StateTwoOld const & /*rhs*/) const;

    StateTwoOld order();

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

#ifndef SWIG
namespace std {

template <>
struct hash<StateOneOld> {
    size_t operator()(const StateOneOld &s) const {
        std::size_t seed = 0;
        utils::hash_combine(seed, s.n);
        utils::hash_combine(seed, s.l);
        utils::hash_combine(seed, s.j);
        utils::hash_combine(seed, s.m);
        return seed;
    }
};

template <>
struct hash<StateTwoOld> {
    size_t operator()(const StateTwoOld &s) const {
        std::size_t seed = 0;
        utils::hash_combine(seed, s.n);
        utils::hash_combine(seed, s.l);
        utils::hash_combine(seed, s.j);
        utils::hash_combine(seed, s.m);
        return seed;
    }
};

} // namespace std
#endif

#endif
