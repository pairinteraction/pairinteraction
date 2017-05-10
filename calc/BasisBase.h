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
    Basis() : energy_min(std::numeric_limits<double>::lowest()), energy_max(std::numeric_limits<double>::max()), range_n({}), range_l({}), range_j({}), range_m({}), hash_restrictions_old(0) {
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
    void restrictJ(double j_min, double j_max) {
        j_min = std::fmax(j_min, 0.5);
        range_j.resize(j_max-j_min+1);
        std::iota(range_j.begin(), range_j.end(), j_min);
        range_j_hash_new = boost::hash_value(range_j);
    }
    void restrictJ(std::vector<double> j) {
        range_j = j;
        range_j_hash_new = boost::hash_value(range_j);
    }
    void restrictM(double m_min, double m_max) {
        range_m.resize(m_max-m_min+1);
        std::iota(range_m.begin(), range_m.end(), m_min);
        range_m_hash_new = boost::hash_value(range_m);
    }
    void restrictM(std::vector<double> m) {
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
        energies = e;
        coefficients = c;
    }

    // TODO void removeUnnecessaryStates(); setThreshold
protected:
    virtual void initialize() = 0;

    double energy_min, energy_max;
    std::vector<int> range_n, range_l;
    std::vector<double> range_j, range_m;

    size_t range_energy_hash, range_n_hash, range_l_hash, range_j_hash, range_m_hash;
    size_t range_energy_hash_new, range_n_hash_new, range_l_hash_new, range_j_hash_new, range_m_hash_new;

    std::vector<double> energies;
    std::vector<T> states;
    eigen_sparse_t coefficients;

private:
    void build() {
        if (energies.empty() && states.empty() && coefficients.size() == 0) {
            initialize();
            range_energy_hash = range_energy_hash_new;
            range_n_hash = range_n_hash_new;
            range_l_hash = range_l_hash_new;
            range_j_hash = range_j_hash_new;
            range_m_hash = range_m_hash_new;
        } else if (!energies.empty() && !states.empty() && coefficients.size() != 0) {
            update();
        }
    }
    void update() {
        if (range_n_hash_new != range_n_hash || range_l_hash_new != range_l_hash ||
                range_j_hash_new != range_j_hash || range_m_hash_new != range_m_hash) {

            std::unordered_set<int> range_set_n;
            if (range_n_hash_new != range_n_hash) range_set_n.insert(range_n.begin(), range_n.end());
            std::unordered_set<int> range_set_l;
            if (range_l_hash_new != range_l_hash) range_set_l.insert(range_l.begin(), range_l.end());
            std::unordered_set<double> range_set_j;
            if (range_j_hash_new != range_j_hash) range_set_j.insert(range_j.begin(), range_j.end());
            std::unordered_set<double> range_set_m;
            if (range_m_hash_new != range_m_hash) range_set_m.insert(range_m.begin(), range_m.end());

            size_t idx_new = 0;

            std::vector<T> states_new;
            states_new.reserve(states.size());
            std::vector<eigen_triplet_double> triplets_transformator;
            triplets_transformator.reserve(states.size());

            for (size_t idx = 0; idx < states.size(); ++idx) {
                const T &state = states[idx];
                if (checkIsQuantumnumberContained(state.n, range_set_n) &&
                        checkIsQuantumnumberContained(state.l, range_set_l) &&
                        checkIsQuantumnumberContained(state.j, range_set_j) &&
                        checkIsQuantumnumberContained(state.m, range_set_m)) {

                    states_new.push_back(state);
                    triplets_transformator.push_back(eigen_triplet_double(idx_new++,idx,1));
                }
            }

            states_new.shrink_to_fit();
            eigen_sparse_t transformator(idx_new,states.size());
            transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

            states = states_new;
            coefficients = transformator*coefficients;

            range_n_hash = range_n_hash_new;
            range_l_hash = range_l_hash_new;
            range_j_hash = range_j_hash_new;
            range_m_hash = range_m_hash_new;
        }

        if (range_energy_hash_new != range_energy_hash ||
                range_n_hash_new != range_n_hash || range_l_hash_new != range_l_hash ||
                range_j_hash_new != range_j_hash || range_m_hash_new != range_m_hash) {
            // Remove cols, and energies based on energy / normalization threshold
            // TODO

            range_energy_hash = range_energy_hash_new;
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


    // SWIG extensions for SciPy interoperability
#ifdef SWIG
    size_t _getCoefficientsNumrows() {
        return getCoefficients().rows();
    }
    size_t _getCoefficientsNumcols() {
        return getCoefficients().cols();
    }
    std::vector<size_t> _getCoefficientsIndptr() {
        std::vector<size_t> data(getCoefficients().outerIndexPtr(), getCoefficients().outerIndexPtr()+getCoefficients().outerSize());
        data.push_back(getCoefficients().nonZeros());
        return data;
    }
    std::vector<size_t> _getCoefficientsIndices() {
        std::vector<size_t> data(getCoefficients().innerIndexPtr(), getCoefficients().innerIndexPtr()+getCoefficients().nonZeros());
        return data;
    }
    std::vector<std::complex<double>> _getCoefficientsData() {
        std::vector<std::complex<double>> data(getCoefficients().valuePtr(), getCoefficients().valuePtr()+getCoefficients().nonZeros());
        return data;
    }
#endif

};


#endif
