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

#ifndef BASIS_H
#define BASIS_H

#include "dtypes.h"
#include "Iter.h"
#include "State.h"
#include "ConfParser.h"

#include <vector>
#include <numeric>
#include <string>
#include <set>
#include <unordered_set>
#include <memory>
#include <string>
#include <fstream>

#include <boost/functional/hash.hpp>

template<class T> class Basis {
public:
    Basis() : energy_min(0), energy_max(0), range_n({}), range_l({}), range_j({}), range_m({}), hash_restrictions_old(0) { // TODO set energy to its limits instead of zero
    }
    void restrictEnergy(double e_min, double e_max) {
        energy_min = e_min;
        energy_max = e_max;
    }
    void restrictN(int n_min, int n_max) {
        n_min = std::max(n_min, 1);
        range_n.resize(n_max-n_min+1);
        std::iota(range_n.begin(), range_n.end(), n_min);
    }
    void restrictN(std::vector<int> n) {
        range_n = n;
    }
    void restrictL(int l_min, int l_max) {
        l_min = std::max(l_min, 0);
        range_l.resize(std::abs(l_max-l_min)+1);
        std::iota(range_l.begin(), range_l.end(), l_min);
    }
    void restrictL(std::vector<int> l) {
        range_l = l;
    }
    void restrictJ(float j_min, float j_max) {
        j_min = std::fmax(j_min, 0.5);
        range_j.resize(j_max-j_min+1);
        std::iota(range_j.begin(), range_j.end(), j_min);
    }
    void restrictJ(std::vector<float> j) {
        range_j = j;
    }
    void restrictM(float m_min, float m_max) {
        range_m.resize(m_max-m_min+1);
        std::iota(range_m.begin(), range_m.end(), m_min);
    }
    void restrictM(std::vector<float> m) {
        range_m = m;
    }
    const std::vector<double>& getEnergies() { // TODO use real_t and make SWIG work with real_t, too
        this->build();
        return energies;
    }
    const std::vector<T>& getStates() {
        this->build();
        return states;
    }
    // TODO void removeUnnecessaryStates();
protected:
    virtual void build() = 0;
    bool checkNewRestrictions() {
        size_t hash_restrictions_new = 0;
        boost::hash_combine(hash_restrictions_new, energy_min);
        boost::hash_combine(hash_restrictions_new, energy_max);
        boost::hash_combine(hash_restrictions_new, range_n);
        boost::hash_combine(hash_restrictions_new, range_l);
        boost::hash_combine(hash_restrictions_new, range_j);
        boost::hash_combine(hash_restrictions_new, range_m);
        if (hash_restrictions_new != hash_restrictions_old) {
            hash_restrictions_old = hash_restrictions_new;
            return true;
        } else {
            return false;
        }
    }
    real_t energy_min, energy_max;
    std::vector<int> range_n, range_l;
    std::vector<float> range_j, range_m;

    std::vector<real_t> energies;
    std::vector<T> states;
    eigen_sparse_t coefficients;

private:
    size_t hash_restrictions_old;
};

class BasisOne : public Basis<StateOne> {
public:
    BasisOne(std::string element) : element(element) {
    }
protected:
    void build() {
        std::vector<int> range_adapted_n, range_adapted_l;
        std::vector<float> range_adapted_j, range_adapted_m;

        if (this->checkNewRestrictions() || this->checkNewElement()) {

            // TODO check whether specified basis is finite

            size_t idx = 0;
            std::vector<eigen_triplet_t> coefficients_triplets; // TODO reserve states, energies, basis_triplets

            if (range_n.empty()) {
                range_adapted_n = std::vector<int>({}); // TODO if empty, calculate the range via the energies
            } else {
                range_adapted_n = range_n;
            }
            for (auto n : range_adapted_n) {

                if (range_l.empty()) {
                    range_adapted_l.resize(n);
                    std::iota(range_adapted_l.begin(), range_adapted_l.end(), 0);
                } else {
                    range_adapted_l = range_l;
                }
                for (auto l : range_adapted_l) {
                    if (l > n-1) continue;

                    if (range_j.empty()) {
                        range_adapted_j = (l == 0) ? std::vector<float>({l+0.5f}) : std::vector<float>({std::fabs(l-0.5f), l+0.5f});
                    } else {
                        range_adapted_j = range_j;
                    }
                    for (auto j : range_adapted_j) {
                        if (std::fabs(j-l) != 0.5) continue;

                        real_t energy = StateOne(element,n,l,j,0.5).getEnergy();
                        if ((energy < energy_min && energy_min != 0) || (energy > energy_max  && energy_max != 0)) continue; // TODO take into account numerical errors

                        if (range_m.empty()) {
                            range_adapted_m.resize(2*j+1);
                            std::iota(range_adapted_m.begin(), range_adapted_m.end(), -j);
                        } else {
                            range_adapted_m = range_m;
                        }
                        for (auto m : range_adapted_m) {
                            if (std::fabs(m) > j) continue;

                            states.push_back(StateOne(element,n,l,j,m));
                            energies.push_back(energy);
                            coefficients_triplets.push_back(eigen_triplet_t(idx,idx,1)); // TODO take into account symmetries

                            ++idx;
                        }
                    }
                }
            }

            coefficients.resize(idx,idx);
            coefficients.setFromTriplets(coefficients_triplets.begin(), coefficients_triplets.end());
            coefficients_triplets.clear();
        }
    }
private:
    bool checkNewElement() {
        return true; // TODO
    }
    std::string element;
};

class BasisTwo : public Basis<StateTwo> {
public:
    BasisTwo(const BasisOne &b1, const BasisOne &b2) : basis1(b1), basis2(b2) {
    }
    BasisOne getFirstBasis() const {
        return basis1;
    }
    BasisOne getSecondBasis() const {
        return basis2;
    }
    void setFirstBasis(const BasisOne &b) {
        basis1 = b;
    }
    void setSecondBasis(const BasisOne &b) {
        basis2 = b;
    }
protected:
    void build() {
        if (this->checkNewRestrictions() || this->checkNewBasisOne()) {
            // Restrict one atom states to the allowed quantum numbers
            std::unordered_set<int> range_set_n(range_n.begin(), range_n.end());
            std::unordered_set<int> range_set_l(range_l.begin(), range_l.end());
            std::unordered_set<float> range_set_j(range_j.begin(), range_j.end());
            std::unordered_set<float> range_set_m(range_m.begin(), range_m.end());

            std::vector<StateOne> states1_restricted; // TODO reserve
            eigen_sparse_t coefficients1_restricted; // TODO fill this variable
            for (const auto &state: basis1.getStates()) {
                if ((range_set_n.empty() || range_set_n.find(state.getN()) != range_set_n.end()) &&
                        (range_set_l.empty() || range_set_l.find(state.getL()) != range_set_l.end()) &&
                        (range_set_j.empty() || range_set_j.find(state.getJ()) != range_set_j.end()) &&
                        (range_set_m.empty() || range_set_m.find(state.getM()) != range_set_m.end())) {
                    states1_restricted.push_back(state);
                }
            }

            std::vector<StateOne> states2_restricted; // TODO reserve
            eigen_sparse_t coefficients2_restricted; // TODO fill this variable
            for (const auto &state: basis2.getStates()) {
                if ((range_set_n.empty() || range_set_n.find(state.getN()) != range_set_n.end()) &&
                        (range_set_l.empty() || range_set_l.find(state.getL()) != range_set_l.end()) &&
                        (range_set_j.empty() || range_set_j.find(state.getJ()) != range_set_j.end()) &&
                        (range_set_m.empty() || range_set_m.find(state.getM()) != range_set_m.end())) {
                    states2_restricted.push_back(state);
                }
            }

            // Combine the one atom states
            // TODO
        }
    }
private:
    bool checkNewBasisOne() {
        return true; // TODO
    }
    BasisOne basis1; // TODO maybe, make const, pass by reference
    BasisOne basis2; // TODO maybe, make const, pass by reference
};

#endif
