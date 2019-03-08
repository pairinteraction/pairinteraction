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

#ifndef BASISNAMES_H
#define BASISNAMES_H

#include "ConfParser.h"
#include "Iter.h"
#include "StateOld.h"
#include "dtypes.h"

#include <fstream>
#include <memory>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

template <class T>
class Basisnames {
public:
    Basisnames() = default;
    void configure(const Configuration &config) {
        conf["deltaNSingle"] << config["deltaNSingle"];
        conf["deltaLSingle"] << config["deltaLSingle"];
        conf["deltaJSingle"] << config["deltaJSingle"];
        conf["deltaMSingle"] << config["deltaMSingle"];
        conf["deltaNSingle"] >> delta_n;
        conf["deltaLSingle"] >> delta_l;
        conf["deltaJSingle"] >> delta_j;
        conf["deltaMSingle"] >> delta_m;
        conf["missingCalc"] << config["missingCalc"];
        conf["missingWhittaker"] << config["missingWhittaker"];
    }
    size_t size() const { return names_.size(); }
    size_t dim() const { return dim_; }
    T &get(size_t idx) { return names_[idx]; }
    const T &get(size_t idx) const { return names_[idx]; }
    void set(size_t i, const T &v) { names_[i] = v; }
    ConstIter<Basisnames, T> begin() const { return ConstIter<Basisnames, T>(this, 0); }
    ConstIter<Basisnames, T> end() const { return ConstIter<Basisnames, T>(this, names_.size()); }
    const Configuration &
    getConf() const { // TODO in Configurable Klasse auslagern, von der geerbt werrden soll
        return conf;
    }

protected:
    int delta_n, delta_l, delta_j, delta_m;
    Configuration conf;
    std::vector<T> names_;
    size_t dim_;
};

class BasisnamesTwo;

class BasisnamesOne : public Basisnames<StateOneOld> {
public:
    /*BasisnamesOne(const Configuration &config, const StateOne &startstate);
    BasisnamesOne(const Configuration &config, const StateTwo &startstate);*/
    static BasisnamesOne fromStates(const std::vector<StateOneOld> &names); // TODO
    static BasisnamesOne fromFirst(const Configuration &config);
    static BasisnamesOne fromFirst(const std::shared_ptr<const BasisnamesTwo> &basis_two);
    static BasisnamesOne fromSecond(const Configuration &config);
    static BasisnamesOne fromSecond(const std::shared_ptr<const BasisnamesTwo> &basis_two);
    static BasisnamesOne fromBoth(const Configuration &config);
    const std::vector<StateOneOld> &initial() const;
    void removeUnnecessaryStates(const std::vector<bool> &is_necessary);
    bool constructedFromFirst();
    void save(const std::string &path);

private:
    BasisnamesOne();
    void build(StateOneOld startstate, const std::string &species);
    void build(StateOneOld startstate, const std::string &species,
               const std::shared_ptr<const BasisnamesTwo> &basis_two, int i);
    void build(StateTwoOld startstate, const std::string &species);
    std::vector<StateOneOld> states_initial;
    bool _constructedFromFirst;
};

class BasisnamesTwo : public Basisnames<StateTwoOld> {
public:
    BasisnamesTwo(const std::shared_ptr<const BasisnamesOne> &basis_one1,
                  const std::shared_ptr<const BasisnamesOne> &basis_one2);
    BasisnamesTwo(const std::shared_ptr<const BasisnamesOne> &basis_one1);
    const StateTwoOld &initial() const;
    void removeUnnecessaryStates(const std::vector<bool> &is_necessary);
    void removeUnnecessaryStatesKeepIdx(const std::vector<bool> &is_necessary);
    void save(const std::string &path);

protected:
    void build(StateTwoOld startstate, std::array<std::string, 2> species,
               const std::shared_ptr<const BasisnamesOne> &basis_one1,
               const std::shared_ptr<const BasisnamesOne> &basis_one2);

private:
    StateTwoOld state_initial;
};

#endif
