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

#ifndef BASISNAMES_H
#define BASISNAMES_H

#include "ConfParser.h"
#include "Iter.h"
#include "State.h"
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

class BasisnamesOne : public Basisnames<StateOne> {
public:
    /*BasisnamesOne(const Configuration &config, const StateOne &startstate);
    BasisnamesOne(const Configuration &config, const StateTwo &startstate);*/
    static BasisnamesOne fromStates(const std::vector<StateOne> &names); // TODO
    static BasisnamesOne fromFirst(const Configuration &config);
    static BasisnamesOne fromFirst(const std::shared_ptr<const BasisnamesTwo> &basis_two);
    static BasisnamesOne fromSecond(const Configuration &config);
    static BasisnamesOne fromSecond(const std::shared_ptr<const BasisnamesTwo> &basis_two);
    static BasisnamesOne fromBoth(const Configuration &config);
    const std::vector<StateOne> &initial() const;
    void removeUnnecessaryStates(const std::vector<bool> &is_necessary);
    bool constructedFromFirst();
    void save(const std::string &path);

private:
    BasisnamesOne();
    void build(StateOne startstate, const std::string &species);
    void build(StateOne startstate, const std::string &species,
               const std::shared_ptr<const BasisnamesTwo> &basis_two, int i);
    void build(StateTwo startstate, const std::string &species);
    std::vector<StateOne> states_initial;
    bool _constructedFromFirst;
};

class BasisnamesTwo : public Basisnames<StateTwo> {
public:
    BasisnamesTwo(const std::shared_ptr<const BasisnamesOne> &basis_one1,
                  const std::shared_ptr<const BasisnamesOne> &basis_one2);
    BasisnamesTwo(const std::shared_ptr<const BasisnamesOne> &basis_one1);
    const StateTwo &initial() const;
    void removeUnnecessaryStates(const std::vector<bool> &is_necessary);
    void removeUnnecessaryStatesKeepIdx(const std::vector<bool> &is_necessary);
    void save(const std::string &path);

protected:
    void build(StateTwo startstate, std::array<std::string, 2> species,
               const std::shared_ptr<const BasisnamesOne> &basis_one1,
               const std::shared_ptr<const BasisnamesOne> &basis_one2);

private:
    StateTwo state_initial;
};

#endif
