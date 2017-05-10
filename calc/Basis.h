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
#include <limits>
#include <unordered_set>
#include <memory>
#include <string>
#include <fstream>

#include <boost/functional/hash.hpp>

template<class T> class Basis {
public:
    Basis() : energy_min(std::numeric_limits<real_t>::lowest()), energy_max(std::numeric_limits<real_t>::max()), range_n({}), range_l({}), range_j({}), range_m({}), hash_restrictions_old(0) {
    }
    void restrictEnergy(double e_min, double e_max) {
        energy_min = e_min;
        energy_max = e_max;
        range_energy_hash_new = boost::hash_value(energy_min);
        boost::hash_combine(range_energy_hash_new, energy_max);
    }
    void restrictN(int n_min, int n_max) {
        n_min = std::max(n_min, 1);
        range_n.resize(n_max-n_min+1);
        std::iota(range_n.begin(), range_n.end(), n_min);
        range_n_hash_new = boost::hash_value(range_n);
    }
    void restrictN(std::vector<int> n) {
        range_n = n;
        range_n_hash_new = boost::hash_value(range_n);
    }
    void restrictL(int l_min, int l_max) {
        l_min = std::max(l_min, 0);
        range_l.resize(std::abs(l_max-l_min)+1);
        std::iota(range_l.begin(), range_l.end(), l_min);
        range_l_hash_new = boost::hash_value(range_l);
    }
    void restrictL(std::vector<int> l) {
        range_l = l;
        range_l_hash_new = boost::hash_value(range_l);
    }
    void restrictJ(float j_min, float j_max) {
        j_min = std::fmax(j_min, 0.5);
        range_j.resize(j_max-j_min+1);
        std::iota(range_j.begin(), range_j.end(), j_min);
        range_j_hash_new = boost::hash_value(range_j);
    }
    void restrictJ(std::vector<float> j) {
        range_j = j;
        range_j_hash_new = boost::hash_value(range_j);
    }
    void restrictM(float m_min, float m_max) {
        range_m.resize(m_max-m_min+1);
        std::iota(range_m.begin(), range_m.end(), m_min);
        range_m_hash_new = boost::hash_value(range_m);
    }
    void restrictM(std::vector<float> m) {
        range_m = m;
        range_m_hash_new = boost::hash_value(range_m);
    }
    const std::vector<double>& getEnergies() { // TODO use real_t and make SWIG work with real_t, too
        this->build();
        return energies;
    }
    const std::vector<T>& getStates() {
        this->build();
        return states;
    }
    const eigen_sparse_t& getCoefficients() {
        this->build();
        return coefficients;
    }
    void updateEnergiesAndCoefficients(const std::vector<real_t> &e, eigen_sparse_t &c) {
        // TODO check whether dimensions fit between e and c and are compatible with states as well
        energies = e;
        coefficients = c;
    }

    // TODO void removeUnnecessaryStates(); setThreshold
protected:
    virtual void initialize() = 0;

    real_t energy_min, energy_max;
    std::vector<int> range_n, range_l;
    std::vector<float> range_j, range_m;

    size_t range_energy_hash, range_n_hash, range_l_hash, range_j_hash, range_m_hash;
    size_t range_energy_hash_new, range_n_hash_new, range_l_hash_new, range_j_hash_new, range_m_hash_new;

    std::vector<real_t> energies;
    std::vector<T> states;
    eigen_sparse_t coefficients;

private:
    void build() {
        // Generate the basis from scratch or update the basis
        if (energies.empty() && states.empty() && coefficients.size() == 0) {
            initialize();
        } else if (!energies.empty() && !states.empty() && coefficients.size() != 0) {
            update();
        }

        // Update hash
        range_energy_hash = range_energy_hash_new;
        range_n_hash = range_n_hash_new;
        range_l_hash = range_l_hash_new;
        range_j_hash = range_j_hash_new;
        range_m_hash = range_m_hash_new;
    }
    void update() {
        if (range_n_hash_new != range_n_hash || range_l_hash_new != range_l_hash ||
                range_j_hash_new != range_j_hash || range_m_hash_new != range_m_hash) {

            // TODO check states.size() == coefficients.innerSize() == coefficients.rows()

            // Transform the vectors of valid quantum numbers into sets
            std::unordered_set<int> range_set_n;
            if (range_n_hash_new != range_n_hash) range_set_n.insert(range_n.begin(), range_n.end());
            std::unordered_set<int> range_set_l;
            if (range_l_hash_new != range_l_hash) range_set_l.insert(range_l.begin(), range_l.end());
            std::unordered_set<float> range_set_j;
            if (range_j_hash_new != range_j_hash) range_set_j.insert(range_j.begin(), range_j.end());
            std::unordered_set<float> range_set_m;
            if (range_m_hash_new != range_m_hash) range_set_m.insert(range_m.begin(), range_m.end());

            // Build transformator and remove states (if the quantum numbers are not allowed)
            std::vector<T> states_new;
            states_new.reserve(states.size());
            std::vector<eigen_triplet_real_t> triplets_transformator;
            triplets_transformator.reserve(states.size());

            size_t idx_new = 0;
            for (size_t idx = 0; idx < states.size(); ++idx) {
                const T &state = states[idx];
                if (checkIsQuantumnumberContained(state.getN(), range_set_n) &&
                        checkIsQuantumnumberContained(state.getL(), range_set_l) &&
                        checkIsQuantumnumberContained(state.getJ(), range_set_j) &&
                        checkIsQuantumnumberContained(state.getM(), range_set_m)) {

                    states_new.push_back(state);
                    triplets_transformator.push_back(eigen_triplet_real_t(idx_new++,idx,1));
                }
            }

            states_new.shrink_to_fit();
            eigen_sparse_real_t transformator(idx_new,states.size());
            transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

            states = states_new;

            // Apply transformator in order to remove rows from the coefficient matrix (i.e. states)
            coefficients = transformator*coefficients;
        }

        if (range_energy_hash_new != range_energy_hash ||
                range_n_hash_new != range_n_hash || range_l_hash_new != range_l_hash ||
                range_j_hash_new != range_j_hash || range_m_hash_new != range_m_hash) {

            // TODO check energies.size() == coefficients.outerSize() == coefficients.cols(),

            // Build transformator and remove energies (if the qenergy is not allowed or the squared norm to small)
            std::vector<double> energies_new;
            energies_new.reserve(energies.size());
            std::vector<eigen_triplet_real_t> triplets_transformator;
            triplets_transformator.reserve(energies.size());

            size_t idx_new = 0;
            for (size_t idx=0; idx<coefficients.outerSize(); ++idx) { // idx = col = num basis vector

                if ((energies[idx] > energy_min || energy_min == std::numeric_limits<real_t>::lowest()) && (energies[idx] < energy_max  || energy_max == std::numeric_limits<real_t>::max())) {
                    real_t sqnorm = 0;
                    for (eigen_iterator_t triple(coefficients,idx); triple; ++triple) {
                        sqnorm += std::pow(std::abs(triple.value()),2);
                    }
                    if (sqnorm > 0.05) { // TODO setThresholdForSQNorm
                        triplets_transformator.push_back(eigen_triplet_real_t(idx,idx_new++,1));
                        energies_new.push_back(energies[idx]);
                    }
                }
            }

            energies_new.shrink_to_fit();
            eigen_sparse_real_t transformator(energies.size(),idx_new);
            transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

            energies = energies_new;

            // Apply transformator in order to remove columns from the coefficient matrix (i.e. basis vectors)
            coefficients = coefficients*transformator;

            // TODO check states.size() == coefficients.innerSize() == coefficients.rows()

            // Build transformator and remove states (if the squared norm is to small)

            std::vector<double> sqnorm_list(states.size(),0);
            for (eigen_idx_t k=0; k<coefficients.outerSize(); ++k) {
                for (eigen_iterator_t triple(coefficients,k); triple; ++triple) {
                    sqnorm_list[triple.row()] += std::pow(std::abs(triple.value()),2);
                }
            }

            std::vector<T> states_new;
            states_new.reserve(states.size());
            triplets_transformator.clear();
            triplets_transformator.reserve(states.size());

            idx_new = 0;
            for (size_t idx = 0; idx < states.size(); ++idx) {
                if (sqnorm_list[idx] > 0.05) {
                    states_new.push_back(states[idx]);
                    triplets_transformator.push_back(eigen_triplet_real_t(idx_new++,idx,1));
                }
            }

            states_new.shrink_to_fit();
            transformator = eigen_sparse_real_t(idx_new,states.size());
            transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

            states = states_new;

            // Apply transformator in order to remove rows from the coefficient matrix (i.e. states)
            coefficients = transformator*coefficients;
        }
    }
    template<class V>
    bool checkIsQuantumnumberContained(V q, std::unordered_set<V> range_q) {
        return range_q.empty() || range_q.find(q) != range_q.end();
    }
    template<class V>
    bool checkIsQuantumnumberContained(std::array<V,2> q, std::unordered_set<V> range_q) {
        return range_q.empty() || (range_q.find(q[0]) != range_q.end() && range_q.find(q[1]) != range_q.end());
    }

    size_t hash_restrictions_old;
};

class BasisOne : public Basis<StateOne> {
public:
    BasisOne(std::string element) : element(element) {
    }
protected:
    void initialize() {
        // TODO check whether specified basis is finite

        size_t idx = 0;
        std::vector<eigen_triplet_t> coefficients_triplets; // TODO reserve states, energies, basis_triplets
        std::vector<int> range_adapted_n, range_adapted_l;
        std::vector<float> range_adapted_j, range_adapted_m;

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
                    if ((energy < energy_min && energy_min != std::numeric_limits<real_t>::lowest()) || (energy > energy_max  && energy_max != std::numeric_limits<real_t>::max())) continue; // TODO take into account numerical errors

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
private:
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
    void initialize() {
        // Restrict one atom states to the allowed quantum numbers
        basis1.restrictN(range_n);
        basis1.restrictL(range_l);
        basis1.restrictJ(range_j);
        basis1.restrictM(range_m);
        basis2.restrictN(range_n);
        basis2.restrictL(range_l);
        basis2.restrictJ(range_j);
        basis2.restrictM(range_m);

        // Combine the one atom states
        // TODO
    }
private:
    bool checkNewBasisOne() {
        return true; // TODO
    }
    BasisOne basis1; // TODO maybe, make const, pass by reference
    BasisOne basis2; // TODO maybe, make const, pass by reference
};

#endif
