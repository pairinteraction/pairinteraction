#ifndef BASIS_BASE_H
#define BASIS_BASE_H

#include "dtypes.h"
#include "State.h"

#include <vector>
#include <numeric>
#include <string>
#include <set>
#include <limits>
#include <unordered_set>
#include <memory>
#include <string>

template<class T> class Basis {
public:
    Basis(std::vector<T> states) : // TODO
        energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()), range_n({}), range_l({}), range_j({}), range_m({}), range_energy_hash(0), range_n_hash(0), range_l_hash(0), range_j_hash(0), range_m_hash(0), states(states) {
    }
    virtual ~Basis() = default;
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
    const std::vector<double>& getEnergies() { // TODO use double and make SWIG work with double, too
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
    void updateEnergiesAndCoefficients(const std::vector<double> &e, eigen_sparse_t &c) {
        // TODO check whether dimensions fit between e and c and are compatible with states as well

        // Set the basis from the input
        energies = e;
        coefficients = c;

        // Reset hash so that the basis is updated accordingly to the energy restrictions in case of access
        range_energy_hash = 0;
    }
    void build() {
        // Generate the basis from scratch or update the basis
        if (energies.empty() && states.empty() && coefficients.size() == 0) { // TODO implement a method that also works if !states.empty()
            initialize();
        } else if (!energies.empty() && !states.empty() && coefficients.size() != 0) {
            update();
        } else {
            // TODO rise error
        }

        // Update hash
        range_energy_hash = range_energy_hash_new;
        range_n_hash = range_n_hash_new;
        range_l_hash = range_l_hash_new;
        range_j_hash = range_j_hash_new;
        range_m_hash = range_m_hash_new;
    }
    void clear() {
        energies.clear();
        states.clear();
        coefficients.resize(0,0);
    }

    // SWIG extensions for SciPy interoperability
    //#ifdef SWIG // ERROR "#ifdef SWIG" does not work for me, here. Definitions seems to be ignored if something is included with %template. In addition, it seems as these methods need to be public (I know, this is not nice ...). I would prefer to write swig related functions within an extra file included by Interfaces.i.cmakein. There, the additional python definitions could be written down, too.
    size_t _getCoefficientsNumrows() {
        return coefficients.rows();
    }
    size_t _getCoefficientsNumcols() {
        return coefficients.cols();
    }
    std::vector<size_t> _getCoefficientsIndptr() {
        std::vector<size_t> data(coefficients.outerIndexPtr(), coefficients.outerIndexPtr()+coefficients.outerSize());
        data.push_back(getCoefficients().nonZeros());
        return data;
    }
    std::vector<size_t> _getCoefficientsIndices() {
        std::vector<size_t> data(coefficients.innerIndexPtr(), coefficients.innerIndexPtr()+coefficients.nonZeros());
        return data;
    }
    std::vector<scalar_t> _getCoefficientsData() {
        std::vector<scalar_t> data(coefficients.valuePtr(), coefficients.valuePtr()+coefficients.nonZeros());
        return data;
    }
    //#endif

    // TODO void removeUnnecessaryStates(); setThreshold
protected:
    Basis() : energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()), range_n({}), range_l({}), range_j({}), range_m({}), range_energy_hash(0), range_n_hash(0), range_l_hash(0), range_j_hash(0), range_m_hash(0) {
    }
    virtual void initialize() = 0;

    double energy_min, energy_max;
    std::vector<int> range_n, range_l;
    std::vector<float> range_j, range_m;

    size_t range_energy_hash, range_n_hash, range_l_hash, range_j_hash, range_m_hash;
    size_t range_energy_hash_new, range_n_hash_new, range_l_hash_new, range_j_hash_new, range_m_hash_new;

    std::vector<double> energies;
    std::vector<T> states;
    eigen_sparse_t coefficients;

private:
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
            std::vector<eigen_triplet_t> triplets_transformator; // TODO use eigen_triplet_int_t
            triplets_transformator.reserve(states.size());

            size_t idx_new = 0;
            for (size_t idx = 0; idx < states.size(); ++idx) {
                const T &state = states[idx];
                if (checkIsQuantumnumberContained(state.n, range_set_n) &&
                        checkIsQuantumnumberContained(state.l, range_set_l) &&
                        checkIsQuantumnumberContained(state.j, range_set_j) &&
                        checkIsQuantumnumberContained(state.m, range_set_m)) {

                    states_new.push_back(state);
                    triplets_transformator.push_back(eigen_triplet_t(idx_new++,idx,1));
                }
            }

            states_new.shrink_to_fit();
            eigen_sparse_t transformator(idx_new,states.size()); // TODO use eigen_sparse_int_t
            transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

            states = states_new;

            // Apply transformator in order to remove rows from the coefficient matrix (i.e. states)
            coefficients = transformator*coefficients;
        }

        if (range_energy_hash_new != range_energy_hash ||
                range_n_hash_new != range_n_hash || range_l_hash_new != range_l_hash ||
                range_j_hash_new != range_j_hash || range_m_hash_new != range_m_hash) {

            // TODO check energies.size() == coefficients.outerSize() == coefficients.cols()

            // Build transformator and remove energies (if the qenergy is not allowed or the squared norm to small)
            std::vector<double> energies_new;
            energies_new.reserve(energies.size());
            std::vector<eigen_triplet_t> triplets_transformator;
            triplets_transformator.reserve(energies.size());

            size_t idx_new = 0;
            for (size_t idx=0; idx<energies.size(); ++idx) { // idx = col = num basis vector

                if ((energies[idx] > energy_min || energy_min == std::numeric_limits<double_t>::lowest()) && (energies[idx] < energy_max  || energy_max == std::numeric_limits<double_t>::max())) {
                    double_t sqnorm = 0;
                    for (eigen_iterator_t triple(coefficients,idx); triple; ++triple) {
                        sqnorm += std::pow(std::abs(triple.value()),2);
                    }
                    if (sqnorm > 0.05) { // TODO setThresholdForSQNorm
                        triplets_transformator.push_back(eigen_triplet_t(idx,idx_new++,1));
                        energies_new.push_back(energies[idx]);
                    }
                }
            }

            energies_new.shrink_to_fit();
            eigen_sparse_t transformator(energies.size(),idx_new);
            transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

            energies = energies_new;

            // Apply transformator in order to remove columns from the coefficient matrix (i.e. basis vectors)
            coefficients = coefficients*transformator;

            // TODO check states.size() == coefficients.innerSize() == coefficients.rows()

            // Build transformator and remove states (if the squared norm is to small)
            std::vector<double> sqnorm_list(states.size(),0);
            for (int k=0; k<coefficients.outerSize(); ++k) {
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
                    triplets_transformator.push_back(eigen_triplet_t(idx_new++,idx,1));
                }
            }

            states_new.shrink_to_fit();
            transformator = eigen_sparse_t(idx_new,states.size());
            transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

            states = states_new;

            // Apply transformator in order to remove rows from the coefficient matrix (i.e. states)
            coefficients = transformator*coefficients;
        }

        // TODO range_quantumnumber after it was used for restrictions, hash is than also not necessary
    }

    template<class V>
    bool checkIsQuantumnumberContained(V q, std::unordered_set<V> range_q) {
        return range_q.empty() || range_q.find(q) != range_q.end();
    }

    template<class V>
    bool checkIsQuantumnumberContained(std::array<V, 2> q, std::unordered_set<V> range_q) {
        return range_q.empty() || (range_q.find(q[0]) != range_q.end() && range_q.find(q[1]) != range_q.end());
    }

};


#endif
